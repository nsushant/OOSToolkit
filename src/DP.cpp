
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
#include "DP.hpp"
#include "Nbody.hpp"



double raan_min(double tstep_size, double a_depot, double a_client, double inclination, double e){

    // this calculates the max possible relative raan rate in rads/sec
    double diff_raan_rate = std::abs(raan_rate(a_depot, e, inclination) - raan_rate(a_client,e,inclination));

    return tstep_size*diff_raan_rate;

}


// Recursion depth protection
double raan_min_res = raan_min(tstep_size, a_depot, a_client, inclination, 0.001);


double calculate_raan_delta_degrees(double raan1_rad, double raan2_rad) {
    double raan1_deg = raan1_rad * 180.0 / M_PI;
    double raan2_deg = raan2_rad * 180.0 / M_PI;
    double diff = std::abs(raan2_deg - raan1_deg);
    return std::min(diff, 360.0 - diff); // Handle wrap-around (359°→1° = 2°)
}



lookupkey make_key(orbital_elements dep_orb, orbital_elements arr_orb, double tof, unsigned int flight_type ){

  lookupkey k ={

        static_cast<int>(std::round(calculate_raan_delta_degrees(dep_orb.RAAN, arr_orb.RAAN) * (raan_min_res * (180.0 / M_PI)) )),
        static_cast<unsigned int>(std::abs(std::round(arr_orb.semi_major_axis - dep_orb.semi_major_axis))),
        flight_type,
        static_cast<unsigned int>(std::abs(std::round((tof) /tstep_size)))

  };

return k;

}



void adjust_trajectories(   schedule_struct& sched, int& blim, double& deltaV_min,
                            std::unordered_map<lookupkey, double, HashKey>& lookup,
                            DataFrame simfile ){

    double deltaV_local = deltaV_min;
    double deltaV_local_add = deltaV_local;
    double deltaV_local_sub = deltaV_local;
    task_block fromblock, toblock, fromblock_sub, toblock_sub;
    double dep_add,dep_sub;
    arma::vec v1sol, v2sol, r1sol, r2sol;
    std::vector<arma::vec> trajs;
    double tof_optimal;

    int b = 1;

    while(b <= blim){

            bool feasibility_condition = false;


            if(b< blim){

                feasibility_condition = ((sched.blocks[b].departure_time <= sched.blocks[b+1].arrival_time) &&
                                         (sched.blocks[b].departure_time <= sched.blocks[b].arrival_constraint));

            }

            else{

                feasibility_condition = ( sched.blocks[b].departure_time <= sched.blocks[b].arrival_constraint );

            }

            if(feasibility_condition == true){

                //positive increment

                deltaV_local_add -= sched.blocks[b].deltaV_arrival;
                deltaV_local_add -= sched.blocks[b+1].deltaV_arrival;

                fromblock = sched.blocks[b];
                toblock = sched.blocks[b+1];

                find_optimal_trajectory(fromblock.satname, toblock.satname, fromblock.departure_time,
                    toblock.arrival_time, v1sol, v2sol, r1sol, r2sol, tof_optimal, trajs, simfile, "edelbaum",
                    false, toblock.deltaV_arrival, lookup);


                deltaV_local_add += toblock.deltaV_arrival;


                dep_add = sched.blocks[b-1].departure_time;

                find_optimal_trajectory_no_iter(sched.blocks[b-1].satname, fromblock.satname, dep_add,
                    fromblock.arrival_time, simfile, fromblock.deltaV_arrival, "edelbaum");

                deltaV_local_add += fromblock.deltaV_arrival;


                //negative increment

                // b - 1 can be at minimum 0 (when b = 1).

                // at b=1 or b-1 =0 we do not have a defined arrival delta V
                deltaV_local_sub -= (b > 1) ? sched.blocks[b - 1].deltaV_arrival: 0.0 ;

                deltaV_local_sub -= sched.blocks[b].deltaV_arrival ;


                fromblock_sub = sched.blocks[b-1];
                toblock_sub = sched.blocks[b];

                find_optimal_trajectory(fromblock_sub.satname, toblock_sub.satname, fromblock_sub.departure_time,
                    toblock_sub.arrival_time, v1sol, v2sol, r1sol, r2sol, tof_optimal, trajs, simfile, "edelbaum",
                    false, toblock_sub.deltaV_arrival);

                deltaV_local_add += toblock_sub.deltaV_arrival;

                //ensuring that the fromblock has a valid previous index
                if(b > 1){

                    dep_sub = sched.blocks[b-2].departure_time;

                    find_optimal_trajectory_no_iter(sched.blocks[b-2].satname, fromblock_sub.satname, dep_sub,
                    fromblock_sub.arrival_time, simfile, fromblock_sub.deltaV_arrival, "edelbaum");
                    deltaV_local_add += fromblock_sub.deltaV_arrival;

                }

                if(std::min(deltaV_local_add,deltaV_local_sub) < deltaV_local){

                    if(deltaV_local_add > deltaV_local_sub){

                        deltaV_local = deltaV_local_sub;
                        sched.blocks[b-1] = fromblock_sub;
                        sched.blocks[b] = toblock_sub;

                    }

                    else{

                        deltaV_local = deltaV_local_add;
                        sched.blocks[b+1] = toblock;
                        sched.blocks[b] = fromblock;

                    }

                }

            }

        b++;

    }

}








void DP (schedule_struct& sched, double& deltav_init, DataFrame simfile){

    int blim = 0;

    double deltaV_min = 1e23;

    std::unordered_map<lookupkey, double, HashKey> lookup;

    while (blim < sched.blocks.size()){

        // update sched block positions to minimium delta V positions.
        // update deltaVmin if improvement is found.
        // if improvement is not found terminate while loop.

        adjust_trajectories(sched, blim, deltaV_min, lookup, simfile);

        blim++;

    }

}
