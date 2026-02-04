#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "LambertSolver.hpp"
#include "Nbody.hpp"



bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2);


void finding_individual_minimas_dynamic_programming(schedule_struct &init_schedule, DataFrame Simfile, double dt);


double calculate_deltav_upto_thisblock(int n, schedule_struct&sched_in, DataFrame simfile, std::map<std::vector<double>, double>&lookup);

int dynamic_program_fixed_tasksize_Tfixed(schedule_struct& sched_init, double dt, int b_reached, int b_lim, double Tmax, DataFrame simfile);


double deltavtotcalc(schedule_struct schedin);

void blocked_time_min(schedule_struct& schedsub, double& time_min);

void branch( schedule_struct& schedule_branch, const int& b, schedule_struct &schedin,
             std::vector<int> visited, double& incumbent, const DataFrame& simfile );

schedule_struct branch_and_bound(schedule_struct schedin, DataFrame simfile, double incumbent);
