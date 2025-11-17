#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "LambertSolver.hpp"


bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2);



schedule_struct finding_individual_minimas_dynamic_programming(schedule_struct init_schedule, DataFrame Simfile, double dt);




