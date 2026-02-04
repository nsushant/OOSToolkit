#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <unordered_set>
#include <functional>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "LambertSolver.hpp"
#include "Nbody.hpp"
#include "LowThrustAnalytical.hpp"

/*
 The outputs of local search algorithms need to be compared with those of exact methods, so as to
 quantify the difference between the solution obtained and the optimal solution (to calculate the optimality gap).

author: Sushanta Nigudkar
date: 11/2025

*/

// Safe access helper functions for bounds checking
inline task_block& safe_block_access(schedule_struct& sched, size_t index, const std::string& context) {
    assert(index < sched.blocks.size() && "Block index out of bounds");
    if (index >= sched.blocks.size()) {
        throw std::out_of_range("Block index " + std::to_string(index) +
                               " out of bounds in " + context);
    }
    return sched.blocks[index];
}

inline const task_block& safe_block_access(const schedule_struct& sched, size_t index, const std::string& context) {
    assert(index < sched.blocks.size() && "Block index out of bounds");
    if (index >= sched.blocks.size()) {
        throw std::out_of_range("Block index " + std::to_string(index) +
                               " out of bounds in " + context);
    }
    return sched.blocks[index];
}

// Recursion depth protection
static const int MAX_RECURSION_DEPTH = 1000;
thread_local int recursion_depth = 0;

// Floating point comparison epsilon
const double EPSILON = 1e-10;

// Transfer delta cache key structure
struct TransferDeltaKey {
    int delta_semi_major_axis;    // Round to nearest meter
    int delta_inclination_deg;    // Round to nearest degree  
    int delta_raan_deg;          // Round to nearest degree
    int transfer_time_bin;        // Round to nearest 1000s
    
    bool operator==(const TransferDeltaKey& other) const {
        return delta_semi_major_axis == other.delta_semi_major_axis &&
               delta_inclination_deg == other.delta_inclination_deg &&
               delta_raan_deg == other.delta_raan_deg &&
               transfer_time_bin == other.transfer_time_bin;
    }
};

// Hash function for unordered_map
struct TransferDeltaKeyHash {
    size_t operator()(const TransferDeltaKey& key) const {
        return std::hash<int>{}(key.delta_semi_major_axis) ^ 
               std::hash<int>{}(key.delta_inclination_deg << 8) ^ 
               std::hash<int>{}(key.delta_raan_deg << 16) ^ 
               std::hash<int>{}(key.transfer_time_bin << 24);
    }
};

// Global transfer delta cache accessible to all recursive branches
static std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash> transfer_delta_cache;

// RAAN wrap-around delta calculation
double calculate_raan_delta_degrees(double raan1_rad, double raan2_rad) {
    double raan1_deg = raan1_rad * 180.0 / M_PI;
    double raan2_deg = raan2_rad * 180.0 / M_PI;
    double diff = std::abs(raan2_deg - raan1_deg);
    return std::min(diff, 360.0 - diff); // Handle wrap-around (359°→1° = 2°)
}

// Cache usage function
double get_or_compute_transfer_deltaV(
    const TransferDeltaKey& key,
    std::function<double()> compute_function
) {
    auto it = transfer_delta_cache.find(key);
    if (it != transfer_delta_cache.end()) {
        return it->second;  // Cache hit
    }
    
    double deltaV = compute_function();  // Cache miss - compute
    transfer_delta_cache[key] = deltaV;
    return deltaV;
}

// Cache clearing function
void clear_transfer_delta_cache() {
    transfer_delta_cache.clear();
}

// Orbital element accessor for enhanced CSV lookup
class OrbitalElementAccessor {
private:
    DataFrame* simfile;
    
    // Cache for frequently accessed (sat_name, time) combinations
    static std::unordered_map<std::string, std::map<double, orbital_elements>> elem_cache;
    
public:
    OrbitalElementAccessor(DataFrame* df) : simfile(df) {}
    
    orbital_elements get_elements_at_time(const std::string& sat_name, double time) {
        // Check local cache first
        auto sat_it = elem_cache.find(sat_name);
        if (sat_it != elem_cache.end()) {
            auto time_it = sat_it->second.find(time);
            if (time_it != sat_it->second.end()) {
                return time_it->second;
            }
        }
        
        // Load from DataFrame (same pattern as current position/velocity access)
        std::vector<std::string> satnames = (*simfile)["name"];
        arma::uvec sat_idxs = find_idxs_of_match(satnames, sat_name);
        
        arma::vec t_sat = simfile->getNumeric("time_s");
        arma::uvec time_idx = arma::find(t_sat == time);
        
        // Extract orbital elements from enhanced CSV
        arma::vec a_sat = simfile->getNumeric("semi_major_axis");
        arma::vec e_sat = simfile->getNumeric("eccentricity");
        arma::vec i_sat = simfile->getNumeric("inclination");
        arma::vec raan_sat = simfile->getNumeric("RAAN");
        arma::vec omega_sat = simfile->getNumeric("arg_periapsis");
        arma::vec nu_sat = simfile->getNumeric("true_anomaly");
        
        orbital_elements elems;
        elems.semi_major_axis = a_sat(time_idx(0));
        elems.eccentricity = e_sat(time_idx(0));
        elems.inclination = i_sat(time_idx(0));
        elems.RAAN = raan_sat(time_idx(0));
        elems.augment_of_periapsis = omega_sat(time_idx(0));
        elems.true_anomaly = nu_sat(time_idx(0));
        
        // Cache result
        elem_cache[sat_name][time] = elems;
        
        return elems;
    }
};

// Static member definition
std::unordered_map<std::string, std::map<double, orbital_elements>> OrbitalElementAccessor::elem_cache;

// The aim is to use dynamic programming to find an optimal schedule.


bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2){

    std::string name1 = block1.satname;
    std::string name2 = block2.satname;


    std::vector<std::string> satnames = simfile["name"];

    arma::uvec sat_idxs_1= find_idxs_of_match(satnames,name1);
    arma::uvec sat_idxs_2 = find_idxs_of_match(satnames,name2);


    arma::vec t_sat = simfile.getNumeric("time_s");

    arma::vec ts_sat1 = t_sat.elem(sat_idxs_1);
    arma::vec ts_sat2 = t_sat.elem(sat_idxs_2);


    arma::uvec sat1_t1 = arma::find(ts_sat1 == t1);
    arma::uvec sat2_t2 = arma::find(ts_sat2 == t2);

    arma::vec x_sat = simfile.getNumeric("x");
    arma::vec y_sat = simfile.getNumeric("y");
    arma::vec z_sat = simfile.getNumeric("z");


    arma::vec tmp1 = x_sat.elem(sat1_t1);
    arma::vec tmp2 = y_sat.elem(sat1_t1);
    arma::vec tmp3 = z_sat.elem(sat1_t1);

    arma::vec3 vec_r1 = {tmp1(0), tmp2(0), tmp3(0)};

    arma::vec tmp4 = x_sat.elem(sat2_t2);
    arma::vec tmp5 = y_sat.elem(sat2_t2);
    arma::vec tmp6 = z_sat.elem(sat2_t2);

    arma::vec3 vec_r2 = {tmp4(0),tmp5(0),tmp6(0)};

    double r1_norm = arma::norm(vec_r1);
    double r2_norm = arma::norm(vec_r2);


    // algorithm 1 D. Izzo "Revisiting Lambert's Problem"
    arma::vec3 vec_c = vec_r2 - vec_r1;

    // vec norms
    double c = arma::norm(vec_c);

    double s = (r1_norm + r2_norm + c) / 2.0;

    double lambda = sqrt(1.0 - c / s);

    double tof = t2 - t1;

    double T = sqrt( 2.0*MU_EARTH / (s*s*s) ) * tof;

    double T_1 = 2.0/3.0 * (1.0 - pow(lambda,3));

    if(T <= T_1){

        return false;
    }

    else{

        return true;

    }

}




void finding_individual_minimas_dynamic_programming(schedule_struct &init_schedule, DataFrame Simfile, double dt){


    // tabulate this table with solutions as they are calculated

    // table[arrival_time][departure_time]
    // lets discretize into dt second intervals
    //
    // if we start with a schedule that lasts 60000 seconds we will get 6000 timesteps
    // out of these 6000 if dt = 100.
    //
    // There will be further reductions due to the time occupied by servicing tasks themselves
    // (times at which we can neither depart nor arrive)
    //
    // Additionally we can check if a solution has a time of flight = T1 and omit these solutions

    // The table will have a maximum of 4 cols (arrival time, departure time, name, delta V)
    // and a maximum of arrival_time_of_last_block/100 rows.


    //int numtimes = std::ceil(init_schedule.blocks[(int)init_schedule.blocks.size() - 1].arrival_time/dt);
    //int numrows = numtimes * numtimes;

    //view_schedule(init_schedule);

    //std::cout << "dt : " << dt << std::endl;
    //std::cout << "numrows : " << numrows << std::endl;

    //arma::mat table_of_sols(numrows,4,arma::fill::zeros);

    // first lets assume the following constraints apply
    // d_i < a_{i+1}
    // a_i < d_i
    // d_i > a_i + s_i
    //
    //
    // If we now consider fixed arrival times b_i to get additional constraints
    //
    // a_i < b_i
    // d_i > b_i + s_i
    //
    // for each solution we would need the block number, deprature time, arrival time, delta V
    //
    // suppose we get
    //
    //
    //Table of sols [a, d, b] = delta V
    //
    std::map<std::vector<double>,double> table_of_sols;

    schedule_struct optimal_schedule;

    double sols = 0;

    for( int b = 1 ; b < init_schedule.blocks.size() ; b++ ){

        task_block block1 = init_schedule.blocks[b-1];
        task_block block2 = init_schedule.blocks[b];

        double d_constraint = block1.departure_time;
        double a_constraint = block2.arrival_time;


        task_block optblock;
        optblock = block2;


        task_block initblock;
        initblock = block1;

        std::vector<double> dopt;
        std::vector<double> aopt;
        std::vector<double> deltaopt;


        for(double a = a_constraint; a > d_constraint; a-=dt){

            for(double d = d_constraint; d < a; d+=dt){

            //std::cout << "d,a : "<< d <<" , "<<a << std::endl;

            if(is_feasible_sol(block1, block2, Simfile, d, a)){


		        std::vector<double> current_prob = {a,d,(double)b};

                bool match_found = (table_of_sols.find(current_prob) != table_of_sols.end());

                double DeltaVsol = 0;

                if (match_found == false){

                    //std::cout << "computing new soln" << "\n";

                    // OLD: find_optimal_trajectory_no_iter(block1.satname, block2.satname, d, a, Simfile, DeltaVsol, "edelbaum" );
                    
                    // NEW: Orbital element based approach with caching
                    OrbitalElementAccessor accessor(&Simfile);
                    orbital_elements dep_orb = accessor.get_elements_at_time(block1.satname, a); // Both at arrival time
                    orbital_elements arr_orb = accessor.get_elements_at_time(block2.satname, a);

                    TransferDeltaKey key = {
                        static_cast<int>(std::round(arr_orb.semi_major_axis - dep_orb.semi_major_axis)),
                        static_cast<int>(std::round((arr_orb.inclination - dep_orb.inclination) * 180.0 / M_PI)),
                        static_cast<int>(std::round(calculate_raan_delta_degrees(dep_orb.RAAN, arr_orb.RAAN))),
                        static_cast<int>(std::round((a - d) / 1000.0))
                    };

DeltaVsol = get_or_compute_transfer_deltaV(key, [&dep_orb, &arr_orb, &a, &d]() {
                        // Use the new orbital element-based function for enhanced performance
                        return calculate_edelbaum_deltaV_orbital_elems(dep_orb, arr_orb, a - d, "passive");
                    });

                    //check if solution exists in lookup table

                    table_of_sols[current_prob] = DeltaVsol;

                    //std::cout <<"wrote solution"<<std::endl;
                    sols+=1;
                }

                else{

                    //std::cout << "trying to read precomputed sol" << "\n";

                    DeltaVsol = table_of_sols[current_prob];

                }


                deltaopt.push_back(DeltaVsol);
                dopt.push_back(d);
                aopt.push_back(a);

            }

          }

        }


        // find time pair that gives minimum delta V
        auto it = std::min_element(deltaopt.begin(), deltaopt.end());
        // Get index
        std::size_t idxmin = std::distance(deltaopt.begin(), it);

        optblock.deltaV_arrival = deltaopt[idxmin];
        optblock.arrival_time = aopt[idxmin];
        initblock.departure_time = dopt[idxmin];


        if (b==1){

            optimal_schedule.blocks.push_back(initblock);

        }


        else{

            optimal_schedule.blocks[ optimal_schedule.blocks.size() - 1 ].departure_time = dopt[idxmin];

        }

        optimal_schedule.blocks.push_back(optblock);

    }


    init_schedule = optimal_schedule;


}





double calculate_deltav_upto_thisblock(int n, schedule_struct&sched_in, DataFrame simfile, std::map<std::vector<double>, double>&lookup){


    if(n==0){

        return 0.0;
    }


    std::vector<double> key = {(double)n,sched_in.blocks[n-1].departure_time,sched_in.blocks[n].arrival_time};


    double deltaV_current = 0.0;

    if (lookup.find(key) != lookup.end()){

        deltaV_current+=lookup[key];

        sched_in.blocks[n].deltaV_arrival = deltaV_current;

    }

    else{

        // if you want to include considered wait-times uncomment the line below

        //deltaV_current += compact_optimal_calc( sched_in.blocks[n-1],sched_in.blocks[n],
            //                                      sched_in[n-1].departure_time, sched_in[n].arrival_time,simfile);


//std::cout <<"n,d,a"<< n <<","<<sched_in.blocks[n-1].departure_time<<","<< sched_in.blocks[n].arrival_time << "\n";

        // OLD: find_optimal_trajectory_no_iter(sched_in.blocks[n-1].satname, sched_in.blocks[n].satname,
        //                                 sched_in.blocks[n-1].departure_time, sched_in.blocks[n].arrival_time,
        //                                 simfile, deltaV_current,"edelbaum");

        // NEW: Orbital element based approach with caching
        OrbitalElementAccessor accessor(&simfile);
        orbital_elements dep_orb = accessor.get_elements_at_time(sched_in.blocks[n-1].satname, sched_in.blocks[n].arrival_time); // Both at arrival time
        orbital_elements arr_orb = accessor.get_elements_at_time(sched_in.blocks[n].satname, sched_in.blocks[n].arrival_time);

        TransferDeltaKey key_transfer = {
            static_cast<int>(std::round(arr_orb.semi_major_axis - dep_orb.semi_major_axis)),
            static_cast<int>(std::round((arr_orb.inclination - dep_orb.inclination) * 180.0 / M_PI)),
            static_cast<int>(std::round(calculate_raan_delta_degrees(dep_orb.RAAN, arr_orb.RAAN))),
            static_cast<int>(std::round((sched_in.blocks[n].arrival_time - sched_in.blocks[n-1].departure_time) / 1000.0))
        };

        deltaV_current = get_or_compute_transfer_deltaV(key_transfer, [&dep_orb, &arr_orb, &sched_in, &n]() {
            // Use the new orbital element-based function for enhanced performance
            return calculate_edelbaum_deltaV_orbital_elems(dep_orb, arr_orb, 
                    sched_in.blocks[n].arrival_time - sched_in.blocks[n-1].departure_time, "passive");
        });

        sched_in.blocks[n].deltaV_arrival = deltaV_current;


        lookup[key] = deltaV_current;

    }

    deltaV_current += calculate_deltav_upto_thisblock(n-1, sched_in,simfile,lookup);

    return deltaV_current;

}



int dynamic_program_fixed_tasksize_Tfixed(schedule_struct& sched_init, double dt,int b_reached,int b_lim,double Tmax,DataFrame simfile,
                                            std::map<std::vector<double>, double>& lookup
){

    // Phase 1: Input validation and bounds checking
    if (sched_init.blocks.empty()) {
        throw std::runtime_error("dynamic_program_fixed_tasksize_Tfixed: Empty schedule not allowed");
    }
    if (b_reached < 0 || b_lim < 0 || b_reached > b_lim) {
        throw std::runtime_error("dynamic_program_fixed_tasksize_Tfixed: Invalid b_reached or b_lim parameters");
    }

    // Phase 1: Recursion depth protection
    if (recursion_depth++ > MAX_RECURSION_DEPTH) {
        recursion_depth--;
        throw std::runtime_error("Maximum recursion depth exceeded in dynamic_program_fixed_tasksize_Tfixed");
    }

    if (b_reached == b_lim){
        recursion_depth--;
        std::cout << " dynamic program terminated" << std::endl;
        return 0;
    }

    //
    // constraints
    // a_i < d_{i-1}
    // a_i = d_i - s_i
    // a_i > 0
    // s_i > 0
    // d_i



    double deltav_so_far = 0.0;

    // Phase 2: Fix off-by-one error - use < instead of <= and bounds checking
    for(int itb = 1 ; itb < b_lim && itb < sched_init.blocks.size(); itb++)
    {
        // Phase 1: Safe array access with bounds checking
        deltav_so_far += safe_block_access(sched_init, itb, "dynamic_program_fixed_tasksize_Tfixed - deltav calculation").deltaV_arrival;
    }


    for(int b = b_reached ; b>=1 ; b--)
    {
        // Phase 1: Bounds checking before array access
        if (b >= sched_init.blocks.size()) {
            std::cerr << "Warning: Skipping invalid b=" << b << " >= blocks.size()=" << sched_init.blocks.size() << std::endl;
            continue;
        }

        std::cout << "b_reached : " << b << std::endl;
        double dep_init = (b > 1) ? safe_block_access(sched_init, b-2, "dynamic_program_fixed_tasksize_Tfixed - dep_init").departure_time : 0.0;

        double arr_init = safe_block_access(sched_init, b, "dynamic_program_fixed_tasksize_Tfixed - arr_init").arrival_time;

        schedule_struct sched_sub = sched_init;

        if (b > 1){

        // Phase 2: Fix floating point comparison with epsilon and bounds checking
        while( safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - while loop").arrival_time > dep_init + EPSILON )
        {

            safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - departure decrement").departure_time -= dt ;

            if( (b-1) > 0 )
            {
                safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - arrival decrement").arrival_time -= dt;

                // Phase 1: Safe array access with bounds checking
                if(safe_block_access(sched_sub, b-2, "dynamic_program_fixed_tasksize_Tfixed - time comparison").departure_time >=
                   safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - time comparison").arrival_time)
                {
                    break;
                }
            }


            // Phase 1: Safe array access for feasibility check
            bool feasibility = is_feasible_sol(safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - feasibility 1"),
                                              safe_block_access(sched_sub, b, "dynamic_program_fixed_tasksize_Tfixed - feasibility 2"),
                                              simfile,
                                              safe_block_access(sched_sub, b-1, "dynamic_program_fixed_tasksize_Tfixed - feasibility 3").departure_time,
                                              safe_block_access(sched_sub, b, "dynamic_program_fixed_tasksize_Tfixed - feasibility 4").arrival_time);

            if(feasibility == false){

                continue;
            }


            double deltav_upto_block = calculate_deltav_upto_thisblock(b_reached, sched_sub, simfile, lookup);


            if (deltav_upto_block < deltav_so_far)
            {
                sched_init = sched_sub;
                deltav_so_far = deltav_upto_block;
            }
        }

        }

        schedule_struct sched_add = sched_init;
        std::cout <<"arr_int"<<arr_init<<std::endl;
        // Phase 2: Fix floating point comparison with epsilon and bounds checking
        while( (safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - second while loop").departure_time < arr_init - EPSILON) &&
               (safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - second while loop constraint").departure_time <
                safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - second while loop constraint").arrival_constraint - EPSILON))
        {

            safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - departure increment").departure_time += dt;
            safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - arrival increment").arrival_time += dt;


            // Phase 1: Safe array access with bounds checking
            if(safe_block_access(sched_add, b, "dynamic_program_fixed_tasksize_Tfixed - time check").arrival_time <=
               safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - time check").departure_time + EPSILON){
                break;
            }




            // Phase 1: Safe array access for feasibility check
            bool feasibility = is_feasible_sol(safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - feasibility 5"),
                                              safe_block_access(sched_add, b, "dynamic_program_fixed_tasksize_Tfixed - feasibility 6"),
                                              simfile,
                                safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - feasibility 7").departure_time,
                                              safe_block_access(sched_add, b-1, "dynamic_program_fixed_tasksize_Tfixed - feasibility 8").arrival_time);

            if(feasibility == false){

                continue;
            }


            double deltav_upto_block = calculate_deltav_upto_thisblock(b_reached, sched_add, simfile, lookup);


            if (deltav_upto_block < deltav_so_far){

                //sched_init.blocks[b-1].departure_time = sched_add.blocks[b-1].departure_time;
                //sched_init.blocks[b-1].arrival_time = sched_add.blocks[b-1].arrival_time;
                sched_init = sched_add;
                deltav_so_far = deltav_upto_block;

            }

        }

    }

    b_reached +=1;

    int result = dynamic_program_fixed_tasksize_Tfixed(sched_init,dt,b_reached,b_lim,Tmax,simfile,lookup);
    recursion_depth--; // Decrement after recursive call
    return result;

}





double deltavtotcalc(schedule_struct schedin){

    double dtot = 0.0;

    for(task_block& block : schedin.blocks){

        dtot+= block.deltaV_arrival;

    }

    return dtot;

}

void blocked_time_min(schedule_struct& schedsub, double& time_min){

    for(task_block block : schedsub.blocks){

        time_min+= block.service_duration;
    }

}


void branch( schedule_struct& schedule_branch, const int& b, schedule_struct &schedin,
             std::vector<int> visited, double& incumbent, const DataFrame& simfile ){
    // Phase 1: Input validation
    if (schedin.blocks.empty()) {
        throw std::runtime_error("branch: Empty schedule not allowed");
    }
    if (b < 0 || static_cast<size_t>(b) >= schedin.blocks.size()) {
        throw std::out_of_range("branch: Invalid block index b = " + std::to_string(b));
    }

    //complete schedule
    if (visited.size() == schedin.blocks.size()){
        schedin = schedule_branch;
        return;
    }

    else{

        // Phase 2: Fix uninitialized variable
        double deltavtot = 0.0;
        schedule_struct bsched;
        bsched = schedule_branch;
        // Phase 1: Safe array access with bounds checking
        bsched.blocks.push_back(safe_block_access(schedin, b, "branch - adding block"));
        std::map<std::vector<double>,double> lookup;

        dynamic_program_fixed_tasksize_Tfixed(  bsched, 100, 0, bsched.blocks.size(),
                                                safe_block_access(bsched, bsched.blocks.size() - 1, "branch - arrival constraint").arrival_constraint,
                                                simfile , lookup );


        deltavtot =  deltavtotcalc(bsched);

        visited.push_back(b);

        std::unordered_set<int> skip(visited.begin(), visited.end());

        double time_min = 0.0;

        blocked_time_min(bsched, time_min);

        double optimistic_estimate = 0.0;



        //calculate bound (optimistic estimate)
        for(int bl = 0; bl < schedin.blocks.size(); bl++){

            if (skip.count(bl)){
                continue;
            }

            // Phase 1: Safe array access with bounds checking
            const task_block& block_bl = safe_block_access(schedin, bl, "branch - optimistic estimate");
            if (time_min > block_bl.arrival_constraint){
                //current branch leaves no time for this block to be
                //feasibly scheduled
                return;
            }

            // Create local non-const copies for the function call
            double time_min_copy = time_min;
            double arrival_constraint_copy = block_bl.arrival_constraint;
            double optimistic_estimate_copy = 0.0; // Function will modify this

            run_exhaustive_search(  safe_block_access(bsched, bsched.blocks.size()-1, "branch - exhaustive search 1").satname,
                                    block_bl.satname,
                                    time_min_copy, arrival_constraint_copy,
                                    optimistic_estimate_copy, "WalkerDelta.csv", "", "edelbaum" );

            optimistic_estimate += optimistic_estimate_copy;


        }




        if( (deltavtot + optimistic_estimate) >= incumbent ){
            // the best case solution for this branch is
            // worse than the incumbent solution
            return;
        }



        for(int blbranch = 0 ; blbranch < schedin.blocks.size() ; blbranch++){

            if (skip.count(blbranch)){
                continue;
            }

            branch( bsched, blbranch, schedin,
                    visited, incumbent, simfile);

        }

    }

}





schedule_struct branch_and_bound(schedule_struct schedin, DataFrame simfile, double incumbent){

  // Phase 1: Input validation
  if (schedin.blocks.empty()) {
    throw std::runtime_error("branch_and_bound: Empty schedule not allowed");
  }
  if (schedin.blocks.size() < 2) {
    throw std::runtime_error("branch_and_bound: Schedule must have at least 2 blocks");
  }

  schedule_struct sol;

  // Safe access with bounds checking
  sol.blocks.push_back( safe_block_access(schedin, 0, "branch_and_bound") );

  std::vector<int> visited;

  visited.push_back(0);
  visited.push_back(schedin.blocks.size() - 1);

  std::unordered_set<int> skip(visited.begin(), visited.end());

  // Phase 2: Fix off-by-one error - ensure minimum size
  if (schedin.blocks.size() < 3) {
    return schedin; // No middle blocks to branch
  }

  for( int b = 1 ; b < schedin.blocks.size() - 1 ; b++ ){

        branch(sol, b, schedin, visited, incumbent, simfile);

  }

  return schedin;

}
