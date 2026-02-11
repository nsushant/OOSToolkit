
#pragma once

#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <functional>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "Trajectory_selection.hpp"



//dynamic program
struct lookupkey{

    // key components
    int delta_raan;
    unsigned int delta_a;
    unsigned int transfer_type; // 1 for sat-sat, 2 for sat-depot, 3 for depot-sat
    unsigned int delta_t;

    // defines the equality operator between two keys
    // the second const just ensures that we do not modify keys when comparing
    // it allows us to check whether two given keys are the same.

    bool operator==(const lookupkey& other) const {

            return delta_a == other.delta_a &&
                   delta_raan == other.delta_raan &&
                   transfer_type == other.transfer_type &&
                   delta_t == other.delta_t;
    }

};



struct HashKey{

    // we need a function that can take in a key and turn it into a single index

    size_t operator()(const lookupkey& key)const{

        return  std::hash<unsigned int>{}(key.delta_a) ^
                ( std::hash<int>{}(key.delta_raan) << 1 ) ^
                (std::hash<unsigned int>{}(key.transfer_type) << 2) ^
                ( std::hash<unsigned int>{}(key.delta_t) << 3 );
    }

};

inline double deltavtotcalc(schedule_struct schedin){

    double dtot = 0.0;

    for(task_block& block : schedin.blocks){

        dtot+= block.deltaV_arrival;

    }

    return dtot;

}


double raan_min(double tstep_size, double a_depot, double a_client, double inclination, double e);


double calculate_raan_delta_degrees(double raan1_rad, double raan2_rad);

lookupkey make_key(orbital_elements dep_orb, orbital_elements arr_orb, double tof, unsigned int flight_type );


void adjust_trajectories(   schedule_struct& sched, int& blim, double& deltaV_min,
                            std::unordered_map<lookupkey, double, HashKey>& lookup,
                            DataFrame simfile );

void DP (schedule_struct& sched, double& deltav_init, DataFrame simfile);
