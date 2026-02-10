#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <unordered_map>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "LambertSolver.hpp"
#include "Nbody.hpp"

// Forward declarations for TransferDeltaKey and its hash function
struct TransferDeltaKey {
    unsigned int delta_semi_major_axis;
    int delta_raan;
    unsigned int transfer_time_bin;
    
    bool operator==(const TransferDeltaKey& other) const {
        return delta_semi_major_axis == other.delta_semi_major_axis &&
               delta_raan == other.delta_raan &&
               transfer_time_bin == other.transfer_time_bin;
    }
};

struct TransferDeltaKeyHash {
    size_t operator()(const TransferDeltaKey& key) const {
        return std::hash<unsigned int>{}(key.delta_semi_major_axis) ^
               (std::hash<int>{}(key.delta_raan) << 1) ^
               (std::hash<unsigned int>{}(key.transfer_time_bin) << 2);
    }
};



bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2);


void finding_individual_minimas_dynamic_programming(schedule_struct &init_schedule, DataFrame Simfile, double dt);


double calculate_deltav_upto_thisblock(int n, schedule_struct&sched_in, DataFrame simfile, std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>&lookup);

int dynamic_program_fixed_tasksize_Tfixed(schedule_struct& sched_init, double dt, int b_reached, int b_lim, double Tmax, DataFrame simfile, std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup);


double deltavtotcalc(schedule_struct schedin);

void blocked_time_min(schedule_struct& schedsub, double& time_min);

void branch( schedule_struct& schedule_branch, const int& b, schedule_struct &schedin,
             std::vector<int> visited, double& incumbent, const DataFrame& simfile );

schedule_struct branch_and_bound(schedule_struct schedin, DataFrame simfile, double incumbent);

// Helper function declarations for refactored DP
bool validate_schedule_constraints(const schedule_struct& schedule, const std::string& context = "");
class ScheduleStateManager;
class DeltaVComparator;
double optimize_block_times_decrement(ScheduleStateManager& state_mgr, int b, int b_reached, 
                                   double dt, double dep_init, DataFrame& simfile,
                                   std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup,
                                   double current_best_deltaV);
double optimize_block_times_increment(ScheduleStateManager& state_mgr, int b, int b_reached,
                                     double dt, double arr_init, DataFrame& simfile,
                                     std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup,
                                     double current_best_deltaV);

// Helper recursive function with shared state manager
int dynamic_program_recursive(schedule_struct& sched_init, ScheduleStateManager& shared_state_mgr, 
                           DeltaVComparator& deltaV_comparator, double dt, int b_reached, int b_lim, double Tmax, DataFrame simfile,
                           std::unordered_map<TransferDeltaKey, double, TransferDeltaKeyHash>& lookup);
