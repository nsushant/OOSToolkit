#include "Local_search.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include <armadillo>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <map>

#include "Local_search.hpp"



/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/

// moves
// use local search to find the optimal_schedule

schedule_struct Tabu_search_opt_schedule_lambert_only(double &init_deltaV, schedule_struct& init_schedule, std::vector<double> dt_move,
                                                        DataFrame simfile, double service_time,std::vector<std::string> move_methods);



void run_tabu_search( DataFrame simfile, std::vector<double> move_size,  
                      std::vector<std::string> moves_to_consider,
                      std::vector<std::string> sat_names_in_schedule ,
                      std::vector<double> t_depart, std::vector<double> t_arrive, 
                      double &deltaV_of_schedule, double service_time);



































