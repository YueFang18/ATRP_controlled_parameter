#include <iomanip>
#include <list>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include "my_rand.h"
#include <time.h>
using namespace std;

double M_total, M_total_square, M_total_cube;
double num_M_chain = 0.0;  // Number of monomers in each polymer chain
int nn = 0;
double num_pr_size = 0.0;

void molecular_weight(double *mol_wt,
                      double &Mn,
                      double &Mw,
                      double *num_M,
                      double *num_Dr,
                      double *num_Pr,
                      double *num_P,
                      vector<int> &Dr,
                      vector<int> &Pr,
                      int *N,
                      double scale_fac,
                      vector<int> poly_M1) {
    // Reset accumulators for molecular weight calculations
    M_total = 0;
    M_total_square = 0;
    *num_M = 0;
    int cnt_poly = 0, cnt_Dr = 0, cnt_Pr = 0;

    // Analyze terminated polymer chains
    for (int i = 0; i < poly_M1.size(); ++i) {
        if (poly_M1[i] != 0) {
            num_M_chain = poly_M1[i];
            *num_M += num_M_chain;
            *num_P += num_M_chain;
            M_total += num_M_chain * (*mol_wt);
            M_total_square += pow(num_M_chain * (*mol_wt), 2);
            nn = poly_M1[i] - 1;
            N[nn]++;
            cnt_poly++;
        }
    }

    // Analyze dormant chains
    for (int i = 0; i < Dr.size(); ++i) {
        if (Dr[i] != 0) {
            num_M_chain = Dr[i];
            *num_M += num_M_chain;
            *num_Dr += num_M_chain;
            M_total += num_M_chain * (*mol_wt);
            M_total_square += pow(num_M_chain * (*mol_wt), 2);
            nn = Dr[i] - 1;
            N[nn]++;
            cnt_Dr++;
        }
    }

    // Analyze active (propagating) chains
    for (int i = 0; i < Pr.size(); ++i) {
        if (Pr[i] != 0) {
            num_M_chain = Pr[i];
            *num_M += num_M_chain;
            *num_Pr += num_M_chain;
            M_total += num_M_chain * (*mol_wt);
            M_total_square += pow(num_M_chain * (*mol_wt), 2);
            num_pr_size = Pr.size();
            nn = Pr[i] - 1;
            N[nn]++;
            cnt_Pr++;
        }
    }

    // Compute molecular weight distributions
    Mn = M_total / (cnt_poly + cnt_Dr + cnt_Pr);           // Number-average molecular weight
    Mw = M_total_square / M_total;                         // Weight-average molecular weight
}
