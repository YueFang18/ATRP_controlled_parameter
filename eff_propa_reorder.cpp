/*
 * main.cpp
 * Hybrid KMC simulation and ODE integration method for ATRP polymerization process
 */

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

// Custom header files
#include "add_schedule.h"                        // Defines ADD_INTERVAL and ADD_TIMES constants
#include "efficient_explicit_sequence_record_number.h" // Updates species counts after reactions
#include "molecular_weight.h"                    // Calculates Mn, Mw
#include "ode.h"                                 // Defines stiff_system for ODE integration

using namespace std;

int main() {
    // Initialize random number generator
    srandom((unsigned)time(NULL));

    // Start timing the simulation
    clock_t start, end;
    start = clock();

    // Set input/output directories
    std::string data_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P2_from_sequence_to_property/effATRP_MMA/src/";
    std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/PhD_project/P2_from_sequence_to_property/effATRP_MMA/output/data20250423/PDI/1e5/20/";

    //-------------------------- Variable Declaration --------------------------//
    // Physical constants
    const double Na = 6.022e23;    // Avogadro's number
    const double R = 0.001986;     // Gas constant (kcal/mol/K)

    // Monomer and polymer properties
    double mol_den;                // Monomer density (g/mL)
    double mol_wt;                 // Monomer molecular weight (g/mol)
    double mol_frac;               // Monomer molar fraction
    double temp;                   // Temperature (K)
    double init_conc[2];           // Initial concentrations of initiator and catalyst (mol/L)
    double monomer_conc;           // Solvent concentration
    double time_limit;             // Maximum simulation time
    double M0;                     // Initial monomer concentration

    // Reaction counters
    int reaction_index;            // Current reaction type index
    int num_chain = -1;            // Chain counter (incremented when new chain created)

    // KMC variables
    double reaction_counter = 0.0; // Cumulative counter for reaction selection
    double tau = 0.0;              // Current reaction time step
    double total_rate = 0.0;       // Sum of all reaction rates
    double reaction = 0.0;         // Random value for reaction selection
    double time = 0.0;             // Current simulation time
    double init_vol = 0.0;         // Initial reaction volume
    double D_time = 0.0;           // Time since last data output
    double d_time = 0.0;           // Time tracker for scaling factor calculation
    double pre_d_time_spe3 = 0.0;  // Time tracker for species 3 (radical) concentration
    double pre_add_spe3 = 0.0;     // Accumulator for radical concentration * time
    double aver_spe3 = 0.0;        // Average radical concentration
    double rand_11, rand_22;       // Random numbers for KMC steps

    // Species concentrations (array index corresponds to species type)
    // [0]: Dr (dormant chains)
    // [1]: C (catalyst)
    // [2]: CX (deactivator)
    // [3]: Pr* (active radical chains)
    // [4]: M (monomer)
    // [5]: P (terminated chains)
    double species[6] = {0};

    // Reaction rates and parameters
    double rate[4] = {0.0};        // Current reaction rates
    double c[4] = {0.0};           // Rate constants

    // Polymerization tracking variables
    double conversion = 0.0;       // Monomer conversion
    double volume = 0.0;           // Reaction volume
    double species_sum = 0.0;      // Total species count (for composition calculation)
    double M_initial = 0.0;        // Initial monomer amount
    double num_monomer = 0.0;      // Monomer consumption counter

    // Molecular weight properties
    double Mn = 0.0;               // Number average molecular weight
    double Mw = 0.0;               // Weight average molecular weight
    double pdi = 0.0;              // Polydispersity index (Mw/Mn)
    double Dpn = 0.0;              // Degree of polymerization

    // Chain counters
    double num_M = 0;              // Monomer units in polymer
    double num_M_sum = 0;          // Total monomer units in polymer
    double num_Dr = 0;             // Dormant chains
    double num_Pr = 0;             // Active chains
    double num_P = 0;              // Dead chains

    // Scaling factor for KMC/ODE hybrid method
    double scale_fac = 1.0;        // Scaling factor for reaction rates
    double radical_ode3 = 0.0;     // Radical concentration from ODE calculation

    // Chain length distribution
    int N[50000] = {0};            // Chain length distribution array

    // Simulation component classes
    poly_chain poly_test;          // Tracks polymer chains
    stiff_system stiff_sys;        // ODE system definition
    stiff_system_jacobi system_jacobi; // Jacobian matrix for ODE system

    // Flags controlling simulation flow
    int flag0 = 0;                 // Tracks first scaling factor calculation
    int flag_add = 0;              // Controls initiator addition

    //-------------------------- Read Input Data --------------------------//
    // Read simulation parameters from input file
    ifstream input(data_dir + "input.txt");
    input >> scale_fac;           // Initial scaling factor
    input >> temp;                // Temperature

    input >> mol_den;             // Monomer density (9.4)
    input >> mol_wt;              // Monomer molecular weight (100.0)
    input >> mol_frac;            // Monomer molar fraction (0.58)

    input >> init_conc[0] >> init_conc[1]; // Initial concentrations [0.023, 0.0115]

    input >> monomer_conc;        // Solvent concentration
    input >> time_limit;          // Maximum simulation time
    input >> M0;                  // Initial monomer concentration

    for (int i = 0; i < 4; i++) {
        input >> c[i];            // Rate constants [43.8, 2.4e7, 1.65e3, 2.82e7]
    }
    input.close();

    //-------------------------- Setup Output Files --------------------------//
    // Open output files for various simulation data
    ofstream conv(out_dir + "conversion.out");            // Monomer conversion
    ofstream molwt(out_dir + "molecularweight.out");      // Molecular weight properties
    ofstream speciesout(out_dir + "species.out");         // Species counts
    ofstream species_con(out_dir + "spe_con.out");        // Species concentrations
    ofstream sf(out_dir + "scaling_factor.out");          // Scaling factor data
    ofstream allchain(out_dir + "all_chain.out");         // Chain length data
    ofstream kmc_log(out_dir + "kmc_add_times.txt");      // Initiator addition times

    /*// Set chain length output thresholds (save polymer data at specific Dpn values)
    double dpnThresholds[] = {10, 20, 30, 40, 50, 60, 80, 100};
    const int N_THRESH = sizeof(dpnThresholds) / sizeof(dpnThresholds[0]);
    bool hasOutput[N_THRESH];
    for (int i = 0; i < N_THRESH; i++) {
        hasOutput[i] = false;      // Initialize all thresholds as not reached
    }*/

    // Set chain length output thresholds (save polymer data at specific conversion values)
    double convThresholds[] = {0.10, 0.15, 0.20, 0.25, 0.30,0.35,0.40,0.45,0.50};
    const int N_THRESH_con = sizeof(convThresholds) / sizeof(convThresholds[0]);
    bool hasOutput_con[N_THRESH_con] = {false};

    // Write headers to output files
    conv << "time" << "\t" << "conversion" << endl;
    molwt << "time" << "\t" << "conversion" << "\t" << "Mn" << "\t" << "\t" << "Mw" << "\t" << "pdi" << "\t" << "Dpn" <<  endl;

    speciesout << "time" << "\t" << "species[0]" << "\t" << "species[1] " << "\t" << "species[2]" << "\t"
               << "species[3] " << "\t" << " species[4]" << "\t" << "species[5] " << endl;
    species_con << "time" << "\t" << "species[0]" << "\t" << "species[1] " << "\t" << "species[2]" << "\t"
               << "species[3] " << "\t" << " species[4]" << "\t" << "species[5] " << endl;

    sf << "time" << "\t"  << "scale_fac" << "\t" << "average kmc radical concentration" << "\t"
       << "ode radical concentration" <<  endl;

    //-------------------------- Initialization --------------------------//
    // Calculate initial volume based on monomer amount
    init_vol = (M0 * mol_frac) / (Na * monomer_conc); // 4.67 = normalized monomer concentration
    volume = init_vol;                        // Set total reaction volume
    cout << "Initial volume: " << volume << endl;

    // Calculate initial molecule numbers for each species
    species[0] = round(init_conc[0] * init_vol * Na + 0.5); // Initiator (Dr)
    species[1] = round(init_conc[1] * init_vol * Na + 0.5); // Catalyst (C)
    species[4] = round(monomer_conc * init_vol * Na + 0.5);         // Monomer (M)

    M_initial = species[4];     // Store initial monomer amount

    // Adjust PDI by adding initiator in multiple steps
    int add_times = ADD_TIMES;     // Number of initiator addition steps (from add_schedule.h)
    int add_species0, add_species1;

    // Calculate molecules to add per step
    add_species0 = round(species[0] / add_times + 0.5);  // Initiator per addition
    add_species1 = round(species[1] / add_times + 0.5);  // Catalyst per addition

    // Start with first portion of initiator and catalyst
    species[0] = add_species0;
    species[1] = add_species1;

    // Initialize dormant chain vector
    vector<int> Dr(species[0], 0);             // Dormant chain lengths (initially all 0)
    std::vector<int> Dr_num(add_species0, 0);  // Dormant chain indices

    // Create map to track all chain lengths (key=chain index, value=chain length)
    std::map<int, int> allchainsinfo;

    // Initialize ODE concentration vector (convert molecule counts to concentrations)
    vector_type r(6);
    for (int i = 0; i < 6; ++i) {
        r[i] = species[i] / (Na * volume);  // Convert to molar concentration
        cout << r[i] << " ";
    }
    cout << endl;

    // Initialize all chains with length 0 and assign chain indices
    for (int i = 0; i < Dr.size(); ++i) {
        allchainsinfo[i] = 0;      // All chains initially have length 0
        Dr_num[i] = i;             // Assign chain indices
    }

    // Output initial species counts
    cout << "Initial species[0]: " << species[0] << endl;
    cout << "Initial species[1]: " << species[1] << endl;
    cout << "Initial species[4]: " << species[4] << endl;
    cout << "Total initiator: " << species[0] * add_times << endl;
    cout << "Total catalyst: " << species[1] * add_times << endl;
    cout << "Total monomer: " << species[4] << endl;

    // Initialize reaction rate constants (propensity functions for KMC)
    // rate[0]: Activation - Dr + C → Pr* + CX
    // rate[1]: Deactivation - Pr* + CX → Dr + C
    // rate[2]: Propagation - Pr* + M → Pr* (longer chain)
    // rate[3]: Termination - Pr* + Pr* → P
    rate[0] = c[0] * species[0] * species[1] / (Na * volume);
    rate[1] = 1 / scale_fac * c[1] / (Na * volume) * species[2] * species[3];
    rate[2] = 1 / scale_fac * c[2] / (Na * volume) * species[3] * species[4];
    rate[3] = 1 / scale_fac / scale_fac * c[3] * 2 / (Na * volume) * species[3] * (species[3] - 1) / 2;

    //-------------------------- KMC Simulation Loop --------------------------//
    // Initialize KMC loop variables
    int count_addtimes = 0;                      // Initiator addition counter
    double next_add_time = ADD_INTERVAL;         // Next initiator addition time

    // Main simulation loop (run until time_limit seconds)
    while (time <= time_limit) {
        //-------------- Adjust PDI by staged initiator addition --------------//
        // Add more initiator at scheduled times to broaden molecular weight distribution
        if (time > next_add_time && count_addtimes < (ADD_TIMES - 1) && flag_add == 0) {
            kmc_log << time << endl;  // Record addition time

            // Add new initiator and catalyst
            species[0] += add_species0;
            species[1] += add_species1;

            count_addtimes++;

            // Add new dormant chains to vector
            Dr.insert(Dr.end(), add_species0, 0);

            // Add new chains to chain information map
            int current_size = allchainsinfo.size();
            for (int i = 0; i < add_species0; ++i) {
                allchainsinfo[current_size + i] = 0;
            }

            cout << "Adding initiator at time: " << time << ", new Dr count: " << species[0] << endl;
            flag_add = 1;
            next_add_time += ADD_INTERVAL;  // Schedule next addition
        }
        flag_add = 0;

        //-------------- Generate random numbers for KMC step --------------//
        std::mt19937_64 rng(std::random_device{}());
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        rand_11 = uni(rng);  // For time step calculation
        rand_22 = uni(rng);  // For reaction selection

        //-------------- KMC Step 1: Calculate time step --------------//
        // Calculate total reaction rate
        total_rate = 0.0;
        for (int i = 0; i < 4; i++) total_rate += rate[i];

        // Calculate time step using KMC formula
        tau = (1 / total_rate) * log(1 / rand_11);

        //-------------- KMC Step 2: Select next reaction --------------//
        reaction = rand_22 * total_rate;
        reaction_counter = rate[0];

        reaction_index = 0;
        // Find which reaction occurs based on random value
        while (reaction_counter < reaction) {
            reaction_index++;
            reaction_counter += rate[reaction_index];
        }

        //-------------- KMC Step 3: Update time and accumulators --------------//
        time += tau;        // Advance simulation time
        D_time += tau;      // Time since last data output
        d_time += tau;      // Time for scaling factor calculation

        // Accumulate species[3]*time product for averaging
        pre_add_spe3 += species[3] * tau;
        pre_d_time_spe3 += tau;

        //-------------- Calculate Scaling Factor (SF) --------------//
        // Initial scaling factor calculation (run once at beginning)
        if (d_time > 0.0001 && flag0 == 0 && pre_d_time_spe3 != 0.0) {
            // Integrate ODE system to get ODE-based radical concentration
            ode(stiff_sys, system_jacobi, r, d_time);
            radical_ode3 = r[3];

            // For first calculation, use initial radical concentration of 1
            aver_spe3 = 1;

            // Scaling factor = KMC concentration / ODE concentration
            scale_fac = aver_spe3 / (Na * volume) / radical_ode3;

            // Output scaling factor data
            sf << time << "  " << scale_fac << " " << aver_spe3 / (Na * volume) << " "
               << radical_ode3 << " " << "1111111" << endl;
            cout << "Initial scale factor: " << scale_fac << endl;

            // Reset accumulators
            d_time = 0.0;
            pre_add_spe3 = 0;
            pre_d_time_spe3 = 0;
            flag0 = 1;  // Mark first SF calculation as completed
        }

        // Regular scaling factor calculation (every 20 seconds)
        if (d_time >= 10.0 && pre_d_time_spe3 != 0.0) {
            // Integrate ODE system
            ode(stiff_sys, system_jacobi, r, d_time);
            radical_ode3 = r[3];

            // Calculate average radical concentration from KMC
            aver_spe3 = pre_add_spe3 / pre_d_time_spe3;

            // Update scaling factor
            scale_fac = (aver_spe3 / (Na * volume)) / radical_ode3;

            // Output scaling factor data
            sf << time << "  " << scale_fac << " " << aver_spe3 / (Na * volume) << " "
               << radical_ode3 << " " << pre_add_spe3 / pre_d_time_spe3 << endl;

            // Reset accumulators
            pre_add_spe3 = 0;
            pre_d_time_spe3 = 0;
            d_time = 0.0;
        }

        //-------------- Periodic Data Output (every 1.0 second) --------------//
        if (D_time >= 1.0) {
            // Output species data (counts and concentrations)
            speciesout << time;
            species_con << time;
            speciesout << " " << species[0] << " " << species[1] << " " << species[2] << " "
                       << species[3] << " " << species[4] << " " << species[5] << " " << endl;
            species_con << " " << species[0] / (Na * volume) << " " << species[1] / (Na * volume) << " "
                        << species[2] / (Na * volume) << " " << species[3] / (Na * volume) << " "
                        << species[4] / (Na * volume) << " " << species[5] / (Na * volume) << " " << endl;

            // Output conversion data
            conversion = num_monomer / M_initial;
            conv << time << "\t" << conversion << endl;

            // Calculate molecular weight properties
            molecular_weight(&mol_wt,
                             Mn,          // Number average molecular weight
                             Mw,          // Weight average molecular weight
                             &num_M,      // Monomer units in polymer
                             &num_Dr,     // Dormant chains
                             &num_Pr,     // Active chains
                             &num_P,      // Dead chains
                             Dr,          // Dormant chain vector
                             Pr,          // Active chain vector
                             N,           // Chain length distribution
                             scale_fac,   // Scaling factor
                             poly_test.chain1); // Dead chain vector

            // Calculate polydispersity index
            pdi = Mw / Mn;

            // Update composition tracking
            species_sum += species[4];
            num_M_sum += num_M;

            // Calculate degree of polymerization
            Dpn = num_M_sum / (poly_test.chain1.size() + Dr.size() + Pr.size());

            // Output molecular weight and composition data
            molwt << time << "\t" << conversion << "\t" << Mn << "\t" << "\t" << Mw << "\t" << "\t"
                  << pdi << "\t" << "\t" << Dpn << endl;

            //-------------- Output Chain Info at Dpn Thresholds --------------//
            /*for (int i = 0; i < N_THRESH; i++) {
                // If Dpn threshold reached and not yet output, write chain data
                if (!hasOutput[i] && (Dpn >= dpnThresholds[i])) {
                    hasOutput[i] = true;  // Mark threshold as output

                    // Generate filename for this threshold
                    int val = (int)dpnThresholds[i];
                    std::string fname = out_dir + "chain_length_" + to_string(val) + ".out";
                    ofstream outChain(fname);

                    // Write all chain lengths to file
                    // First dormant chains (Dr)
                    for (size_t k = 0; k < Dr.size(); ++k) {
                        outChain << Dr[k] << " ";
                    }

                    // Then dead chains (poly_test.chain1)
                    for (size_t k = 0; k < poly_test.chain1.size(); ++k) {
                        outChain << poly_test.chain1[k] << " ";
                    }

                    // Finally active chains (Pr)
                    for (size_t k = 0; k < Pr.size(); ++k) {
                        outChain << Pr[k] << " ";
                    }
                    outChain << endl;
                    outChain.close();
                }
            }*/


            //-------------- Output Chain Info at Dpn Thresholds --------------//
            for (int i = 0; i < N_THRESH_con; i++) {
                    if (!hasOutput_con[i] && conversion >= convThresholds[i]) {
                        hasOutput_con[i] = true;
                        cout << "Outputting at conversion: " << convThresholds[i] << endl; // 调试输出
                        // 现有输出代码保持不变
                        std::string fname = out_dir + "all_chain_conv_" + to_string(int(convThresholds[i] * 100)) + ".out";
                        ofstream outChain(fname);
                        if (!outChain.is_open()) {
                            cerr << "Failed to open file: " << fname << endl;
                            continue;
                        }
                        for (size_t k = 0; k < Dr.size(); ++k) {
                            outChain << Dr[k] << " ";
                        }
                        for (size_t k = 0; k < poly_test.chain1.size(); ++k) {
                            outChain << poly_test.chain1[k] << " ";
                        }
                        for (size_t k = 0; k < Pr.size(); ++k) {
                            outChain << Pr[k] << " ";
                        }
                        outChain << endl;
                        outChain.close();
                    }
                }

            // Reset counters for next data collection period
            num_M = 0;
            num_M_sum = 0.0;
            D_time = 0;
        }

        //-------------- KMC Step 4: Execute selected reaction --------------//
        // Update species counts based on selected reaction
        efficient_explicit_sequence_record_number(
                reaction_index,  // Type of reaction to execute
                &num_monomer,    // Monomer consumption counter
                species,         // Species counts
                num_chain,       // Chain counter
                Dr,              // Dormant chain vector
                Dr_num,          // Dormant chain indices
                Pr,              // Active chain vector
                allchainsinfo,   // Chain length mapping
                allchains,       // All chains vector
                poly_test        // Polymer chain data
        );

        //-------------- KMC Step 5: Update reaction rates --------------//
        // Recalculate propensity functions based on updated species counts
        // Activation: Dr + C → Pr* + CX
        rate[0] = c[0] * species[0] * species[1] / (Na * volume);
        rate[1] = 1 / scale_fac * c[1] / (Na * volume) * species[2] * species[3];
        rate[2] = 1 / scale_fac * c[2] / (Na * volume) * species[3] * species[4];
        rate[3] = 1 / scale_fac / scale_fac * c[3] * 2 / (Na * volume) * species[3] * (species[3] - 1) / 2;
    }

    //-------------------------- Finish Simulation --------------------------//
    // Output final chain information and statistics
    cout << "Dr.size(): " << Dr.size() << endl;         // Report number of dormant chains
    cout << "poly_M1.size(): " << poly_test.chain1.size() << endl;  // Report number of dead chains
    cout << "Pr.size(): " << Pr.size() << endl;         // Report number of active chains

    // Write all final chain lengths to output file
    // Format: [dormant chains] [dead chains] [active chains]
    for (size_t i = 0; i < Dr.size(); ++i) {
        allchain << Dr[i] << " ";                     // Dormant chain lengths
    }

    for (size_t i = 0; i < poly_test.chain1.size(); ++i) {
        allchain << poly_test.chain1[i] << " ";       // Dead chain lengths
    }

    for (size_t i = 0; i < Pr.size(); ++i) {
        allchain << Pr[i] << " ";                     // Active chain lengths
    }
    allchain << endl;

    // Calculate and output total run time
    end = clock();
    cout << "Running time = " << (double)(end-start)/CLOCKS_PER_SEC << "s" << endl;

    return 0;
}