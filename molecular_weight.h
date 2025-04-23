#include<iomanip>
#include<list>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include "my_rand.h"
#include<time.h>
using namespace std;

double M_total, M_total_square, M_total_cube;
double num_M_chain[2] = {0.0};  //number of monomers in each chain
//int N[319]={0};
int nn=0;
double num_pr_size=0.0;

void molecular_weight(double *mol_wt,
                      double &Mn,
                      double &Mw,
                      double &Mz,
                      double *num_M,
                      double *num_Dr,
                      double *num_Pr,
                      double *num_P,
                      vector<int> &Dr,
                      vector<int> &Pr,
                      int *N,
                      double scale_fac,
                      vector<int> poly_M1){
  M_total = 0;
  M_total_square = 0;
  M_total_cube = 0;    // 新增变量，用于计算Mz
  num_M[0]=0;
  int cnt_poly = 0, cnt_Dr = 0, cnt_Pr = 0;


  for (int i = 0; i < poly_M1.size(); i++) {
      if (poly_M1[i] !=0){
          num_M_chain[0] = poly_M1[i];
          num_M[0] += num_M_chain[0];
          num_P[0] += num_M_chain[0];
          M_total += num_M_chain[0] * mol_wt[0] ;
          M_total_square += pow((num_M_chain[0] * mol_wt[0]), 2);
          M_total_cube += pow((num_M_chain[0] * mol_wt[0]), 3);
          nn = poly_M1[i];
          nn--;
          N[nn]++;
          cnt_poly++;
      }
  }

    for (int i = 0; i < Dr.size(); i++) {
        if (Dr[i] !=0){
            num_M_chain[0] = Dr[i];
            num_M[0] += num_M_chain[0];
            num_Dr[0] += num_M_chain[0];
            M_total += num_M_chain[0] * mol_wt[0] ;
            M_total_square += pow((num_M_chain[0] * mol_wt[0]), 2);
            M_total_cube += pow((num_M_chain[0] * mol_wt[0]), 3);
            nn = Dr[i];
            nn--;
            N[nn]++;
            cnt_Dr++;
        }

    }

    for (int i = 0; i < Pr.size(); i++) {
        if (Pr[i]!=0){
            num_M_chain[0] = Pr[i];
            num_M[0] += num_M_chain[0];
            num_Pr[0] += num_M_chain[0];
            //M_total += num_M_chain[0] * mol_wt[0] /scale_fac;
            M_total += num_M_chain[0] * mol_wt[0] ;
            //M_total_square += pow((num_M_chain[0] * mol_wt[0])/scale_fac, 2);
            M_total_square += pow((num_M_chain[0] * mol_wt[0]), 2);
            M_total_cube += pow((num_M_chain[0] * mol_wt[0]), 3);
            num_pr_size=Pr.size();
            nn = Pr[i];
            nn--;
            N[nn]++;
            cnt_Pr++;
        }


    }
    //Mn = M_total / (poly_M1.size()+ Dr.size()+ (Pr.size()/scale_fac));
    //Mn = M_total / (poly_M1.size()+ Dr.size()+ (Pr.size()));
    Mn = M_total /(cnt_poly + cnt_Dr + cnt_Pr);
    Mw = M_total_square / M_total;
    Mz = M_total_cube / M_total_square;    // 新增代码行，用于计算Mz
}

