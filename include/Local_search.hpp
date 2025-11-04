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

    double arrival_time; 
    double departure_time;
    std::string satname; 
    double deltaV_arrival; 
    double service_duration;
    double demand_deadline;

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
                                                        double dt_move, DataFrame simfile, double service_time,
                                                        std::vector<std::string> move_methods );




schedule_struct local_search_opt_schedule(double init_deltaV, schedule_struct init_schedule, double dt_move, DataFrame simfile);


