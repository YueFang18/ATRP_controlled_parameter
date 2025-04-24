#include <cstdio>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <time.h>

#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include "add_schedule.h"

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

// Type aliases for the ODE state vector and Jacobian matrix
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

// Reaction rate constants: k0, k1, k2, k3
// Units consistent with mol/L and seconds
static const double k[4] = {43.8, 2.4e7, 1.65e+03, 2.82e7};

// Definition of the stiff reaction system: computes dx/dt = f(x)
struct stiff_system {
    void operator()(const vector_type &x, vector_type &dxdt, double /* t */) {
        // Activation/deactivation equilibrium between dormant (x[0], x[1]) and active species (x[2], x[3])
        dxdt[0] = -k[0] * x[0] * x[1] + k[1] * x[2] * x[3];
        dxdt[1] = -k[0] * x[0] * x[1] + k[1] * x[2] * x[3];
        dxdt[2] =  k[0] * x[0] * x[1] - k[1] * x[2] * x[3];
        // Termination of active radicals (second-order recombination)
        dxdt[3] =  k[0] * x[0] * x[1] - k[1] * x[2] * x[3] - 2 * k[3] * x[3] * x[3];
        // Consumption of monomer by propagation
        dxdt[4] = -k[2] * x[4] * x[3];
        // Formation of terminated products
        dxdt[5] = 2 * k[3] * x[3] * x[3];
    }
};

// Definition of the Jacobian for the stiff system: J = df/dx
struct stiff_system_jacobi {
    void operator()(const vector_type &x, matrix_type &J, const double &/* t */, vector_type &dfdt) {
        // Partial derivatives for species x[0]
        J(0, 0) = -k[0] * x[1];
        J(0, 1) = -k[0] * x[0];
        J(0, 2) =  k[1] * x[3];
        J(0, 3) =  k[1] * x[2];
        // No direct dependence on monomer or product species
        J(0, 4) = 0.0;
        J(0, 5) = 0.0;

        // Rows for x[1] are identical to x[0]
        J(1, 0) = J(0, 0);
        J(1, 1) = J(0, 1);
        J(1, 2) = J(0, 2);
        J(1, 3) = J(0, 3);
        J(1, 4) = 0.0;
        J(1, 5) = 0.0;

        // Partial derivatives for x[2]
        J(2, 0) =  k[0] * x[1];
        J(2, 1) =  k[0] * x[0];
        J(2, 2) = -k[1] * x[3];
        J(2, 3) = -k[1] * x[2];
        J(2, 4) = 0.0;
        J(2, 5) = 0.0;

        // Partial derivatives for x[3] including termination term
        J(3, 0) =  k[0] * x[1];
        J(3, 1) =  k[0] * x[0];
        J(3, 2) = -k[1] * x[3];
        J(3, 3) = -k[1] * x[2] - 4 * k[3] * x[3];
        J(3, 4) = 0.0;
        J(3, 5) = 0.0;

        // Partial derivatives for monomer consumption x[4]
        J(4, 0) = 0.0;
        J(4, 1) = 0.0;
        J(4, 2) = 0.0;
        J(4, 3) = -k[2] * x[4];
        J(4, 4) = -k[2] * x[3];
        J(4, 5) = 0.0;

        // Partial derivatives for product formation x[5]
        J(5, 0) = 0.0;
        J(5, 1) = 0.0;
        J(5, 2) = 0.0;
        J(5, 3) = 4 * k[3] * x[3];
        J(5, 4) = 0.0;
        J(5, 5) = 0.0;

        // No explicit time dependence
        for (size_t i = 0; i < dfdt.size(); ++i) {
            dfdt[i] = 0.0;
        }
    }
};

/**
 * ODE driver synchronized with KMC: integrates stiff system over dt
 * Inputs:
 *   r   - State vector (6 species concentrations in mol/L)
 *   dt  - Time step from KMC (seconds)
 * Output:
 *   r   - Updated state at time t + dt
 **/
void ode(stiff_system, stiff_system_jacobi, vector_type &r, double dt) {
    // Global simulation time and injection count
    static double T_GLOBAL   = 0.0;    // Absolute time (s)
    static int    N_INJECTED = 0;      // Number of injections applied

    // Concentration increment per injection (mol/L), matching KMC initial values
    const double ADD_DR = 0.023  / ADD_TIMES;   // Dormant radical feed
    const double ADD_C  = 0.0115 / ADD_TIMES;   // Monomer feed

    // End time for this integration segment
    const double t_end = T_GLOBAL + dt;
    vector_type x = r;  // Work copy of state

    // Integrate in sub-intervals split by feed events
    while (T_GLOBAL < t_end) {
        // Determine next feed event or segment end
        const double t_next_inj  = (N_INJECTED + 1) * ADD_INTERVAL;
        const double t_stage_end = min(t_next_inj, t_end);

        // Perform stiff integration with Rosenbrock4
        integrate_const(
                make_dense_output<rosenbrock4<double>>(1e-6, 1e-6),
                make_pair(stiff_system(), stiff_system_jacobi()),
                x,
                T_GLOBAL,
                t_stage_end,
                (t_stage_end - T_GLOBAL) / 100.0  // Use 100 substeps per segment
        );

        // Advance global clock
        T_GLOBAL = t_stage_end;

        // At feed events, add concentration increments
        if (T_GLOBAL >= t_next_inj - 1e-6 && N_INJECTED < ADD_TIMES - 1) {
            x[0] += ADD_DR;  // Add dormant radical
            x[1] += ADD_C;   // Add monomer
            ++N_INJECTED;
            // Log injection time
            ofstream ode_log("/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P2_from_sequence_to_property/effATRP_MMA/output/data20250423/ode_add_times.txt", ios::app);
            ode_log << T_GLOBAL << "\n";
            ode_log.close();
        } else if (T_GLOBAL >= t_next_inj - 1e-6) {
            // Final segment: count injection but no concentration addition
            ++N_INJECTED;
        }
    }

    // Return updated state to caller
    r = x;
}
