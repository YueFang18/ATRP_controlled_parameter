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
#include <map>

using namespace std;

struct chain {
    vector<int> M1;
    vector<int> mon_type;
    vector<int> mon_num;
};

struct poly_chain {
    poly_chain(){};
    vector<int> chain1;  //polymer_chian
    vector<int> chain2; //Dr
    vector<int> num_M1, num_M2 ;
} poly;


typedef vector<chain> chain_pool;
chain_pool allchains;                   //completed segments and dead chains
chain new_chain;
int r1, r2, rmax;                       //random numbers for chain selection

vector<int>::iterator p_chain;
int c_tab1, c_tab2;               //2chain labels

std::vector<int> Pr;
std::vector<int> Pr_num;

// 初始化所有链的长度为0，并初始化Dr_num



void efficient_explicit_sequence_record_number(int reaction_index,
                                               double *num_monomer,
                                               double *species,
                                               int &num_chain,
                                               vector<int> &Dr,
                                               vector<int> &Dr_num,
                                               vector<int> &Pr,
                                               map<int, int> &allchainsinfo,
                                               chain_pool &allchains,
                                               poly_chain &poly) {


    switch (reaction_index) {
        /************** Initiator Decomposition **************/
        case 0:              //Dr + C------> Pr* +CX

            species[0] -= 1;
            species[1] -= 1;
            species[2] += 1;
            species[3] += 1;

            r1 = my_rand(Dr.size());
            Pr.push_back(Dr[r1]);
 //           Pr_num.push_back(Dr_num[r1]);  // 在Pr_num中添加这个链的序号

            p_chain = Dr.begin() + r1;
            Dr.erase(p_chain);
//            Dr_num.erase(Dr_num.begin() + r1);  // 在Dr_num中删除这个链的序号
            break;


        case 1:             //Pr*+ CX----> Dr + C

            species[2] -= 1;
            species[3] -= 1;
            species[0] += 1;
            species[1] += 1;

            r1 = my_rand(Pr.size());
            Dr.push_back(Pr[r1]);
//            Dr_num.push_back(Pr_num[r1]);  // 在Dr_num中添加这个链的序号

            p_chain = Pr.begin() + r1;
            Pr.erase(p_chain);
//            Pr_num.erase(Pr_num.begin() + r1);  // 在Pr_num中删除这个链的序号
            break;


            /**************** Chain Propagation*************/
        case 2 :             //Pr* + M -----> Pr+1*

            species[3] -= 1;
            species[3] += 1;
            species[4] -= 1;
            num_monomer[0] += 1;
            r1 = my_rand(Pr.size());
            Pr[r1]++;
 //           allchainsinfo[Pr_num[r1]]++;  // 对应的链在allchains中也增加1
            break;

            /****************Disproportination of type Pr* --> dead end polymer (1 to 2 radical transfer)*************/
        case 3 :       //Pr* + Ps* -----> Pr + Ps

            species[3] -= 2;
            species[5] += 2;

            do {
                    r1 = my_rand(Pr.size());
                    r2 = my_rand(Pr.size());

            } while (r1 == r2);

            c_tab1 = Pr[r1];
            c_tab2 = Pr[r2];

            poly.chain1.push_back(c_tab1);
            poly.chain1.push_back(c_tab2);


            if (r1 < r2) {
                rmax = r2;
                r2 = r1;
                r1 = rmax;
            }

            p_chain = Pr.begin() + r1;
            Pr.erase(p_chain);
//            Pr_num.erase(Pr_num.begin() + r1);  // 在Pr_num中删除这个链的序号

            p_chain = Pr.begin() + r2;
            Pr.erase(p_chain);
//            Pr_num.erase(Pr_num.begin() + r2);  // 在Pr_num中删除这个链的序号
            break;

    }
}
