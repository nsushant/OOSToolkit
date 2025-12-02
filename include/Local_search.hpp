#pragma once 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"

/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local search. 

*/

//structs 
struct task_block{

    double arrival_time= 0.0; 
    double departure_time = 0.0;
    std::string satname = "unassigned"; 
    double deltaV_arrival =0.0; 
    double service_duration =0.0;
    double demand_deadline =0.0;

};

struct schedule_struct{


    std::vector<task_block> blocks;


};

//moves 
void move_dt(task_block &tb, double dt);

void move_dt2(std::vector<task_block> &blocks, int b_index, double dt);
void move_add_arrival(std::vector<task_block> &blocks, int b_index, double dt);
void move_add_departure(std::vector<task_block> &blocks, int b_index, double dt);
void move_sub_departure(std::vector<task_block> &blocks, int b_index, double dt);
void move_sub_arrival(std::vector<task_block> &blocks, int b_index, double dt);

void swap_slots(std::vector<task_block> &blocks, int b_index, double dt,DataFrame simfile);

void move_wrapper(std::vector<task_block> &blocks, int b_index, double dt, std::string method,DataFrame simfile);


void view_schedule(schedule_struct schedule_to_print); 


void init_dep_arrival_times_random(std::vector<double> &departure_times, std::vector<double> &arrival_times,double period, double service_time, int num_sats);
void init_dep_arrival_times_strict_timespan(std::vector<double> &departure_times, std::vector<double> &arrival_times, double final_time, double service_time, int num_sats);
std::vector<std::string> init_satname_array(std::string service_satname, std::vector<std::string> client_satnames, bool allow_revisits, double num_sats);


schedule_struct create_schedule(double &deltaV_of_schedule_init, 
                                std::vector<double> arrival_times,
                                std::vector<double> departure_times, 
                                std::vector<std::string> satnames,
                                DataFrame simfile);



schedule_struct create_schedule_lambert_only(   double &deltaV_of_schedule_init, 
                                                std::vector<double> arrival_times,
                                                std::vector<double> departure_times, 
                                                std::vector<std::string> satnames,
                                                DataFrame simfile, 
                                                double service_time=0.0, 
                                                std::vector<double> deadlines={});



schedule_struct local_search_opt_schedule_lambert_only( double &init_deltaV, schedule_struct init_schedule, 
                                                        std::vector<double> dt_move, DataFrame simfile, double service_time,
                                                        std::vector<std::string> move_methods );




schedule_struct local_search_opt_schedule(double init_deltaV, schedule_struct init_schedule, double dt_move, DataFrame simfile);



schedule_struct run_local_search( DataFrame simfile, std::vector<double> move_size,  
                      std::vector<std::string> moves_to_consider,
                      std::vector<std::string> sat_names_in_schedule ,
                      std::vector<double> t_depart, std::vector<double> t_arrive, 
                      double &deltaV_of_schedule, double service_time);
