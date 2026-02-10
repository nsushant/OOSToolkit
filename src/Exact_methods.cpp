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
// 1. this function accepts non const schedule
inline task_block& safe_block_access(schedule_struct& sched, size_t index, const std::string& context) {
    assert(index < sched.blocks.size() && "Block index out of bounds");
    if (index >= sched.blocks.size()) {
        throw std::out_of_range("Block index " + std::to_string(index) +
                               " out of bounds in " + context);
    }
    return sched.blocks[index];
}

// 2. this function accepts non modifiable const schedules
inline const task_block& safe_block_access(const schedule_struct& sched, size_t index, const std::string& context) {
    assert(index < sched.blocks.size() && "Block index out of bounds");
    if (index >= sched.blocks.size()) {
        throw std::out_of_range("Block index " + std::to_string(index) +
                               " out of bounds in " + context);
    }
    return sched.blocks[index];
}


// calculate binsize for dynamic program

double raan_min(double tstep_size, double a_depot, double a_client, double inclination, double e){

    // this calculates the max possible relative raan rate in rads/sec
    double diff_raan_rate = std::abs(raan_rate(a_depot, e, inclination) - raan_rate(a_client,e,inclination));

    return tstep_size*diff_raan_rate;

}



// Recursion depth protection
double raan_min_res = raan_min(tstep_size, a_depot, a_client, inclination, 0.001);
static const int MAX_RECURSION_DEPTH = 1000;
thread_local int recursion_depth = 0;



// Transfer delta cache key structure
struct TransferDeltaKey {
    unsigned int delta_semi_major_axis;    // Round to nearest meter, always positive
    int delta_raan;                        // RAAN in raan_min_res
    unsigned int transfer_time_bin;        // Round to nearest tstep_size, always positive

    // defines the equality operator between two keys
    // the second const just ensures that we do not modify keys when comparing
    bool operator==(const TransferDeltaKey& other) const {
        return delta_semi_major_axis == other.delta_semi_major_axis &&
               delta_raan == other.delta_raan &&
               transfer_time_bin == other.transfer_time_bin;
    }

};

// Hash function for unordered_map
struct TransferDeltaKeyHash {

    // gives you the index of the key in a hash table
    // you would get the undex by doing TransferDeltaKeyHash h , index = h(key)

    size_t operator()(const TransferDeltaKey& key) const {

        return  std::hash<unsigned int>{}(key.delta_semi_major_axis) ^
               (std::hash<int>{}(key.delta_raan) << 1) ^
               (std::hash<unsigned int>{}(key.transfer_time_bin) << 2);
    }

};

// Global transfer delta cache accessible to all recursive branches
// <Raw_keys, DeltaV_type,hashfunction>
static std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash> transfer_delta_cache;

// RAAN wrap-around delta calculation
double calculate_raan_delta_degrees(double raan1_rad, double raan2_rad) {
    double raan1_deg = raan1_rad * 180.0 / M_PI;
    double raan2_deg = raan2_rad * 180.0 / M_PI;
    double diff = std::abs(raan2_deg - raan1_deg);
    return std::min(diff, 360.0 - diff); // Handle wrap-around (359°→1° = 2°)
}

// Cache usage function
double get_or_compute_transfer_deltaV(const TransferDeltaKey& key, std::function<double()> compute_function)
{
    // see if this transfer has been calculated before
    auto it = transfer_delta_cache.find(key);
    if (it != transfer_delta_cache.end()) {
        return it->second;  // Cache hit
    }

    //else calculated deltaV
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

// Forward declarations for helper functions
bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2);
double calculate_deltav_upto_thisblock(int n, schedule_struct& sched_in, DataFrame simfile,
                                     std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup);

// Constraint validation function for departure_time <= arrival_constraint
bool validate_schedule_constraints(const schedule_struct& schedule, const std::string& context = "") {
    for (size_t i = 0; i < schedule.blocks.size(); ++i) {
        const task_block& block = safe_block_access(schedule, i, "validate_schedule_constraints");

        // Check departure_time <= arrival_constraint
        if (block.departure_time >= block.arrival_constraint ) {
            if (!context.empty()) {
                std::cout << "Constraint violation in " << context << ": Block " << i
                         << " departure_time (" << block.departure_time
                         << ") > arrival_constraint (" << block.arrival_constraint << ")" << std::endl;
            }
            return false;
        }

        // Check arrival_time > departure_time (basic time consistency)
        if (block.arrival_time <= block.departure_time - EPSILON) {
            if (!context.empty()) {
                std::cout << "Time inconsistency in " << context << ": Block " << i
                         << " arrival_time (" << block.arrival_time
                         << ") <= departure_time (" << block.departure_time << ")" << std::endl;
            }
            return false;
        }
    }

    // Check consecutive block timing constraints
    for (size_t i = 1; i < schedule.blocks.size(); ++i) {
        const task_block& prev_block = safe_block_access(schedule, i-1, "validate_schedule_constraints - prev");
        const task_block& curr_block = safe_block_access(schedule, i, "validate_schedule_constraints - curr");

        // arrival_i <= departure_{i-1} should not happen (maintain sequence)
        if (curr_block.arrival_time > prev_block.departure_time ) {
            if (!context.empty()) {
                std::cout << "Sequence violation in " << context << ": Block " << i
                         << " arrival_time (" << curr_block.arrival_time
                         << ") < Block " << (i-1) << " departure_time (" << prev_block.departure_time << ")" << std::endl;
            }
            return false;
        }
    }

    return true;
}

// Schedule state manager for efficient in-place modifications with rollback capability
class ScheduleStateManager {
private:
    schedule_struct& schedule;
    std::vector<std::pair<size_t, task_block>> saved_blocks; // (index, original_block)

public:
    explicit ScheduleStateManager(schedule_struct& sched) : schedule(sched) {}

    // Save block state before modification
    void save_block(size_t index) {
        if (index < schedule.blocks.size()) {
            saved_blocks.emplace_back(index, schedule.blocks[index]);
        }
    }

    // Modify block times
    void modify_block_times(size_t index, double departure_delta, double arrival_delta) {
        if (index < schedule.blocks.size()) {
            save_block(index);
            schedule.blocks[index].departure_time += departure_delta;
            schedule.blocks[index].arrival_time += arrival_delta;
        }
    }

    // Set absolute block times
    void set_block_times(size_t index, double new_departure, double new_arrival) {
        if (index < schedule.blocks.size()) {
            save_block(index);
            schedule.blocks[index].departure_time = new_departure;
            schedule.blocks[index].arrival_time = new_arrival;
        }
    }

    // Rollback all changes
    void rollback() {
        for (const auto& [index, original_block] : saved_blocks) {
            if (index < schedule.blocks.size()) {
                schedule.blocks[index] = original_block;
            }
        }
        saved_blocks.clear();
    }

    // Commit changes (clear saved state)
    void commit() {
        saved_blocks.clear();
    }

    // Get reference to schedule (for read operations)
    const schedule_struct& get_schedule() const { return schedule; }
    schedule_struct& get_schedule() { return schedule; }
};

// Optimize block times by decrementing (moving earlier)
double optimize_block_times_decrement(ScheduleStateManager& state_mgr, int b, int b_reached,
                                   double dt, double dep_init, DataFrame& simfile,
                                   std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup,
                                   double current_best_deltaV) {

    schedule_struct& schedule = state_mgr.get_schedule();

    while (safe_block_access(schedule, b-1, "optimize_block_times_decrement").arrival_time > dep_init + EPSILON) {

        // Apply time decrement
        state_mgr.modify_block_times(b-1, -dt, -dt);

        // Check previous block constraint if applicable
        if ((b-1) > 0) {
            const task_block& prev_block = safe_block_access(schedule, b-2, "optimize_block_times_decrement - prev");
            const task_block& curr_block = safe_block_access(schedule, b-1, "optimize_block_times_decrement - curr");

            if (prev_block.departure_time <= curr_block.arrival_time - EPSILON) {
                break; // Violate sequence constraint
            }
        }

        // Validate all constraints
        if (!validate_schedule_constraints(schedule, "decrement_optimization")) {
            continue;
        }

        // Check feasibility
        bool feasibility = is_feasible_sol(
            safe_block_access(schedule, b-1, "optimize_block_times_decrement - feasibility 1"),
            safe_block_access(schedule, b, "optimize_block_times_decrement - feasibility 2"),
            simfile,
            safe_block_access(schedule, b-1, "optimize_block_times_decrement - feasibility 3").departure_time,
            safe_block_access(schedule, b, "optimize_block_times_decrement - feasibility 4").arrival_time
        );

        if (!feasibility) {
            continue;
        }

        // Calculate deltaV for this configuration
        double new_deltaV = calculate_deltav_upto_thisblock(b_reached, schedule, simfile, lookup);

        // Update if better solution found
        if (new_deltaV < current_best_deltaV) {
            current_best_deltaV = new_deltaV;
            state_mgr.commit(); // Keep the changes
        } else {
            state_mgr.rollback(); // Discard the changes
        }
    }

    return current_best_deltaV;
}

// Optimize block times by incrementing (moving later)
double optimize_block_times_increment(ScheduleStateManager& state_mgr, int b, int b_reached,
                                    double dt, double arr_init, DataFrame& simfile,
                                    std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup,
                                    double current_best_deltaV) {

    schedule_struct& schedule = state_mgr.get_schedule();
    const task_block& current_block = safe_block_access(schedule, b-1, "optimize_block_times_increment");

    while ( (current_block.departure_time < arr_init - EPSILON) &&
            (current_block.departure_time <= current_block.arrival_constraint)) {

        // Apply time increment
        state_mgr.modify_block_times(b-1, dt, dt);

        // Check next block constraint
        const task_block& next_block = safe_block_access(schedule, b, "optimize_block_times_increment - next");
        if (next_block.arrival_time >= current_block.departure_time + EPSILON) {
            break; // Violate sequence constraint
        }

        // Validate all constraints
        if (!validate_schedule_constraints(schedule, "increment_optimization")) {
            continue;
        }

        // Check feasibility
        bool feasibility = is_feasible_sol(
            safe_block_access(schedule, b-1, "optimize_block_times_increment - feasibility 1"),
            safe_block_access(schedule, b, "optimize_block_times_increment - feasibility 2"),
            simfile,
            safe_block_access(schedule, b-1, "optimize_block_times_increment - feasibility 3").departure_time,
            safe_block_access(schedule, b, "optimize_block_times_increment - feasibility 4").arrival_time
        );

        if (!feasibility) {
            continue;
        }

        // Calculate deltaV for this configuration
        double new_deltaV = calculate_deltav_upto_thisblock(b_reached, schedule, simfile, lookup);

        // Update if better solution found
        if (new_deltaV < current_best_deltaV) {
            current_best_deltaV = new_deltaV;
            state_mgr.commit(); // Keep the changes
        } else {
            state_mgr.rollback(); // Discard the changes
        }
    }

    return current_best_deltaV;
}

// Efficient deltaV comparison with memoization
class DeltaVComparator {
private:
    std::unordered_map<size_t, double> deltaV_cache; // key: schedule hash, value: deltaV
    DataFrame& simfile;
    std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& transfer_lookup;

    // Simple hash function for schedule state (based on timing)
    size_t hash_schedule(const schedule_struct& schedule) const {
        size_t hash = 0;
        for (const auto& block : schedule.blocks) {
            hash ^= std::hash<double>{}(block.departure_time) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            hash ^= std::hash<double>{}(block.arrival_time) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }

public:
    DeltaVComparator(DataFrame& sim, std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup)
        : simfile(sim), transfer_lookup(lookup) {}

    // Calculate deltaV with caching
    double calculate_deltaV(int b_reached, schedule_struct& schedule) {
        size_t schedule_hash = hash_schedule(schedule);

        auto it = deltaV_cache.find(schedule_hash);
        if (it != deltaV_cache.end()) {
            return it->second; // Cache hit
        }

        // Cache miss - calculate and store
        double deltaV = calculate_deltav_upto_thisblock(b_reached, schedule, simfile, transfer_lookup);
        deltaV_cache[schedule_hash] = deltaV;
        return deltaV;
    }

    // Compare and update if better
    bool is_better_solution(int b_reached, schedule_struct& schedule, double current_best_deltaV) {
        double new_deltaV = calculate_deltaV(b_reached, schedule);
        return new_deltaV < current_best_deltaV;
    }

    // Clear cache when needed
    void clear_cache() {
        deltaV_cache.clear();
    }
};

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





double calculate_deltav_upto_thisblock(int n, schedule_struct&sched_in, DataFrame simfile, std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>&lookup){


    if(n==0){

        return 0.0;
    }


    // NEW: Orbital element based approach with caching
    OrbitalElementAccessor accessor(&simfile);
    orbital_elements dep_orb = accessor.get_elements_at_time(sched_in.blocks[n-1].satname, sched_in.blocks[n].arrival_time); // Both at arrival time
    orbital_elements arr_orb = accessor.get_elements_at_time(sched_in.blocks[n].satname, sched_in.blocks[n].arrival_time);

    TransferDeltaKey key = {
        static_cast<unsigned int>(std::abs(std::round(arr_orb.semi_major_axis - dep_orb.semi_major_axis))),
        static_cast<int>(std::round(calculate_raan_delta_degrees(dep_orb.RAAN, arr_orb.RAAN) * (raan_min_res * (180.0 / M_PI)) )),
        static_cast<unsigned int>(std::abs(std::round((sched_in.blocks[n].arrival_time - sched_in.blocks[n-1].departure_time) /tstep_size)))
    };

    double deltaV_current = 0.0;

    if (lookup.find(key) != lookup.end()){

        deltaV_current+=lookup[key];

        sched_in.blocks[n].deltaV_arrival = deltaV_current;

    }

    else{


        deltaV_current = get_or_compute_transfer_deltaV(key, [&dep_orb, &arr_orb, &sched_in, &n]() {
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



int dynamic_program_fixed_tasksize_Tfixed(schedule_struct& sched_init, double dt, int b_reached, int b_lim, double Tmax, DataFrame simfile,
                                            std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup) {

    // Input validation
    if (sched_init.blocks.empty()) {
        throw std::runtime_error("dynamic_program_fixed_tasksize_Tfixed: Empty schedule not allowed");
    }
    if (b_reached < 0 || b_lim < 0 || b_reached > b_lim) {
        throw std::runtime_error("dynamic_program_fixed_tasksize_Tfixed: Invalid b_reached or b_lim parameters");
    }

    // Recursion depth protection
    if (recursion_depth++ > MAX_RECURSION_DEPTH) {
        recursion_depth--;
        throw std::runtime_error("Maximum recursion depth exceeded in dynamic_program_fixed_tasksize_Tfixed");
    }

    // Base case
    if (b_reached == b_lim) {
        recursion_depth--;
        std::cout << "Dynamic program terminated successfully" << std::endl;
        return 0;
    }

    // Initialize helper classes
    ScheduleStateManager state_mgr(sched_init);
    DeltaVComparator deltaV_comparator(simfile, lookup);

    // Calculate current best deltaV
    double current_best_deltaV = 0.0;

    for (int itb = 1; itb < b_lim && itb < sched_init.blocks.size(); itb++) {
        current_best_deltaV += safe_block_access(sched_init, itb, "dynamic_program_fixed_tasksize_Tfixed - deltaV init").deltaV_arrival;
    }

    // Optimize each block in reverse order
    for (int b = b_reached; b >= 1; b--) {
        if (b >= sched_init.blocks.size()) {
            std::cerr << "Warning: Skipping invalid b=" << b << " >= blocks.size()=" << sched_init.blocks.size() << std::endl;
            continue;
        }

        std::cout << "Optimizing block: " << b << std::endl;

        // Calculate reference times
        double dep_init = (b > 1) ? safe_block_access(sched_init, b-2, "dynamic_program_fixed_tasksize_Tfixed - dep_init").departure_time : 0.0;
        double arr_init = safe_block_access(sched_init, b, "dynamic_program_fixed_tasksize_Tfixed - arr_init").arrival_time;

        // Phase 1: Decrement optimization (move earlier)
        if (b > 1) {
            current_best_deltaV = optimize_block_times_decrement(state_mgr, b, b_reached, dt, dep_init, simfile, lookup, current_best_deltaV);
        }

        // Phase 2: Increment optimization (move later)
        current_best_deltaV = optimize_block_times_increment(state_mgr, b, b_reached, dt, arr_init, simfile, lookup, current_best_deltaV);
    }

    std::cout<<"best delta V:"<<current_best_deltaV<<"\n";
    // Recursive call for next block
    b_reached += 1;
    int result = dynamic_program_fixed_tasksize_Tfixed(sched_init, dt, b_reached, b_lim, Tmax, simfile, lookup);
    recursion_depth--;
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
        std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash> lookup;

        dynamic_program_fixed_tasksize_Tfixed(  bsched, 1000.0, 0, bsched.blocks.size(),
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
                                    optimistic_estimate_copy, "../data/WalkerDelta.csv", "", "edelbaum" );

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
