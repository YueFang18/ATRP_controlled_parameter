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


/* -----------------------------------------------------------------
 * ODE 驱动 —— 与 KMC 同步的刚性体系积分
 * 传入 :  r  当前时刻的 6 维浓度 (mol·L⁻¹)
 *         dt 本次 KMC 步长 (s)
 * 返回 :  r  更新到 t + dt 的浓度
 * ----------------------------------------------------------------- */
void ode( stiff_system ,
          stiff_system_jacobi ,
          vector_type &r ,
          double dt )
{
    /* ----------- 全局时间 & 已加料次数 ----------- */
    static double T_GLOBAL   = 0.0;    // 绝对时间 (s)
    static int    N_INJECTED = 0;      // 已注料次数

    /* ----------- 每次投料增加的浓度 (mol·L⁻¹) ----------- */
    const double ADD_DR = 0.023  / ADD_TIMES;   // DR  (与 KMC 初值保持一致)
    const double ADD_C  = 0.0115 / ADD_TIMES;   // C

    const double t_end = T_GLOBAL + dt;         // 本次积分的终点

    vector_type x = r;                          // 工作副本

    while( T_GLOBAL < t_end )
    {
        /* --- 本阶段需要积分到的时间 --- */
        const double t_next_inj  = (N_INJECTED + 1) * ADD_INTERVAL;
        const double t_stage_end = std::min( t_next_inj , t_end );

        /* --- Rosenbrock4 刚性积分 --- */
        integrate_const(
                make_dense_output< rosenbrock4<double> >( 1e-6 , 1e-6 ),
                make_pair( stiff_system() , stiff_system_jacobi() ),
                x ,
                T_GLOBAL ,
                t_stage_end ,
                ( t_stage_end - T_GLOBAL ) / 100.0               // 步长 = 1/10 阶段
        );

        T_GLOBAL = t_stage_end;              // 推进全局时钟

        /* --- 跨过投料节点：补加浓度 --- */
        if( T_GLOBAL >= t_next_inj - 1e-6 &&
            N_INJECTED < ADD_TIMES - 1 )
        {
            x[0] += ADD_DR;                  // DR
            x[1] += ADD_C;                   // C
            ++N_INJECTED;
            std::ofstream ode_log("/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P2_from_sequence_to_property/effATRP_MMA/output/data20250413/PDI/1e5/4/ode_add_times.txt", std::ios::app);
            ode_log << T_GLOBAL << "\n";
            ode_log.close();
    }
        else if( T_GLOBAL >= t_next_inj - 1e-6 )
        {
            ++N_INJECTED;                    // 最后一段仅记录次数
        }
    }

    r = x;                                   // 把积分后的状态传回主程序
}



