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
#include "Local_search.hpp"

/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/

double prob_calculation(double old_state_energy,double new_state_energy,double T);


schedule_struct simulated_annealing_lambert(double &init_deltaV, schedule_struct init_schedule, std::vector<double> dt_move,
                                            DataFrame simfile, double service_time,std::vector<std::string> move_methods, double cooling_param, int maxiter, double T_init);

 
void run_simulated_annealing( DataFrame simfile, std::vector<double> move_size,  
                              std::vector<std::string> moves_to_consider,
                              std::vector<std::string> sat_names_in_schedule ,
                              std::vector<double> t_depart, std::vector<double> t_arrive, 
                              double &deltaV_of_schedule, double service_time, double cooling_param, int maxiter, double T_init);


































