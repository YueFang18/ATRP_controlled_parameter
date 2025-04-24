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
    vector<int> chain1;  // Dead polymer chains
    vector<int> chain2;  // Dormant chains
    vector<int> num_M1;  // Monomer counts
} poly;

typedef vector<chain> chain_pool;
chain_pool allchains;                   // Completed chain segments and dead chains
chain new_chain;
int r1, r2, rmax;                       // Random numbers for chain selection

vector<int>::iterator p_chain;
int c_tab1, c_tab2;                     // Chain labels

std::vector<int> Pr;                    // Active chains

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
        /************** Initiator Decomposition (Activation) **************/
        case 0:              // Dr + C → Pr* + CX
            species[0] -= 1;
            species[1] -= 1;
            species[2] += 1;
            species[3] += 1;

            r1 = my_rand(Dr.size());
            Pr.push_back(Dr[r1]);

            p_chain = Dr.begin() + r1;
            Dr.erase(p_chain);
            break;

        case 1:             // Pr* + CX → Dr + C (Deactivation)
            species[2] -= 1;
            species[3] -= 1;
            species[0] += 1;
            species[1] += 1;

            r1 = my_rand(Pr.size());
            Dr.push_back(Pr[r1]);

            p_chain = Pr.begin() + r1;
            Pr.erase(p_chain);
            break;

            /**************** Chain Propagation *************/
        case 2:             // Pr* + M → Pr+1*
            species[3] -= 1;
            species[3] += 1;
            species[4] -= 1;
            *num_monomer += 1;

            r1 = my_rand(Pr.size());
            Pr[r1]++;
            break;

            /**************** Termination Reaction (Radical Coupling) *************/
        case 3:             // Pr* + Ps* → Pr + Ps
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

            // Ensure deletion starts from the higher index to avoid invalidation
            if (r1 < r2) {
                rmax = r2;
                r2 = r1;
                r1 = rmax;
            }

            p_chain = Pr.begin() + r1;
            Pr.erase(p_chain);

            p_chain = Pr.begin() + r2;
            Pr.erase(p_chain);
            break;
    }
}