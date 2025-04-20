/*
 * ODE part
 */
#include <cstdio>
#include<iomanip>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include<time.h>

#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include "add_schedule.h"

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

//[ stiff_system_definition
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;


double k[4] =
        {43.8, 2.4e7, 1.65e+03, 2.82e7}; //2.4e+07,2.82e+07

struct stiff_system {
    void operator()(const vector_type &x, vector_type &dxdt, double /* t */) {
        dxdt[0] = -k[0] * x[0] * x[1] + k[1] * x[2] * x[3];

        dxdt[1] = -k[0] * x[0] * x[1] + k[1] * x[2] * x[3];

        dxdt[2] =  k[0] * x[0] * x[1] - k[1] * x[2] * x[3];

        dxdt[3] =  k[0] * x[0] * x[1] - k[1] * x[2] * x[3] - 2 * k[3] * x[3] * x[3] ;

        dxdt[4] = -k[2] * x[4] *x[3];

        dxdt[5] = 2 * k[3] * x[3] * x[3];

    }
};

struct stiff_system_jacobi {
    void operator()(const vector_type &x, matrix_type &J, const double & /* t */ , vector_type &dfdt) {
        J(0, 0) = -k[0] * x[1];
        J(0, 1) = -k[0] * x[0];
        J(0, 2) = k[1] *  x[3];
        J(0, 3) = k[1] *  x[2];
        J(0, 4) = 0.0;
        J(0, 5) = 0.0;

        J(1, 0) = -k[0] * x[1];
        J(1, 1) = -k[0] * x[0];
        J(1, 2) = k[1] *  x[3];
        J(1, 3) = k[1] *  x[2];
        J(1, 4) = 0.0;
        J(1, 5) = 0.0;

        J(2, 0) = k[0] * x[1];
        J(2, 1) = k[0] * x[0];
        J(2, 2) = - k[1] *  x[3];
        J(2, 3) = - k[1] *  x[2];
        J(2, 4) = 0;
        J(2, 5) = 0;

        J(3, 0) =  k[0] * x[1];
        J(3, 1) =  k[0] * x[0];
        J(3, 2) = - k[1] *  x[3];
        J(3, 3) = - k[1] *  x[2]- 4 * k[3] * x[3];
        J(3, 4) = 0;
        J(3, 5) = 0;

        J(4, 0) = 0.0;
        J(4, 1) = 0;
        J(4, 2) = 0;
        J(4, 3) = -k[2] * x[4];
        J(4, 4) = -k[2] * x[3];
        J(4, 5) = 0.0;

        J(5, 0) = 0.0;
        J(5, 1) = 0;
        J(5, 2) = 0;
        J(5, 3) = 4 * k[3] * x[3];
        J(5, 4) = 0;
        J(5, 5) = 0;

        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
        dfdt[4] = 0.0;
        dfdt[5] = 0.0;

    }
};
//]



/*void ode(stiff_system,stiff_system_jacobi,vector_type &r,vector_type &Rad, double end_t)
{


    std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P2_from_sequence_to_property/effATRP_MMA/output/data20230918/PDI/";

    //double start_t=0.0;
    //double dt_t=5.0;

    vector_type x(6);
    x=r;
    //cout<<r[0]<<endl;

    ofstream t1(out_dir+"ode_data.out");
    size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
                                           make_pair( stiff_system() , stiff_system_jacobi() ) ,
                                           x , 0.0 , end_t ,end_t/10.0,
                                           t1 << phoenix::arg_names::arg2 << " "  << phoenix::arg_names::arg1[0] << " "<< phoenix::arg_names::arg1[1]<< " " << phoenix::arg_names::arg1[2] << " " << phoenix::arg_names::arg1[3] << " " << phoenix::arg_names::arg1[4] << " " << phoenix::arg_names::arg1[5] <<"\n" );

    r=x;
}*/


void ode(stiff_system, stiff_system_jacobi,
         vector_type &r, vector_type & /* Rad */,
         double dt)                       // dt = KMC这一步走了多长时间
{
    /* --------- 静态量只初始化一次 --------- */
    static double INIT_DR = 0.0, INIT_C = 0.0;
    static bool   INIT_SET = false;
    static int    GLOBAL_ADDED = 0;
    static double T_GLOBAL = 0.0;          // <<< 绝对时间基准

    if (!INIT_SET) {
        INIT_DR  = r[0];
        INIT_C   = r[1];
        INIT_SET = true;
    }

    const double add_conc_0 = INIT_DR / ADD_TIMES;
    const double add_conc_1 = INIT_C  / ADD_TIMES;

    double t_end = T_GLOBAL + dt;          // 本次需要积分到的绝对时间
    vector_type x = r;

    while (T_GLOBAL < t_end) {

        // 下一个“该加料”的绝对时刻
        double t_next_inj = (GLOBAL_ADDED + 1) * ADD_INTERVAL;

        // 本阶段积分到 min(下一次加料, 本次时间终点)
        double t_stage_end = std::min(t_next_inj, t_end);

        integrate_const(
                make_dense_output< rosenbrock4<double> >(1e-6, 1e-6),
                make_pair(stiff_system(), stiff_system_jacobi()),
                x,
                T_GLOBAL,                   // 起点
                t_stage_end,                // 终点
                (t_stage_end - T_GLOBAL)/10.0
        );

        T_GLOBAL = t_stage_end;

        /* ----------- 只有真正跨过1000 s倍数才加料 ----------- */
        // 建议的末尾判断（更安全）
        if( T_GLOBAL >= t_next_inj - 1e-12 &&
            GLOBAL_ADDED < ADD_TIMES - 1 )
        {
            x[0] += add_conc_0;
            x[1] += add_conc_1;
            ++GLOBAL_ADDED;
        }

        else if (T_GLOBAL >= t_next_inj - 1e-12) {
            ++GLOBAL_ADDED;   // 记录最后一次，但不再加料
        }
    }

    r = x;   // 把最新状态传回主程序
}


