#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>

#include "Local_search.hpp"
#include "data_access_lib.hpp"

/*
 The outputs of local search methods find only local minima, to calidate thier
accuracy their solutions have to be compared to those reached by exact methods
like branch and bound or branch and cut. This file implements an exact algorithm
for this purpose.

author: Sushanta Nigudkar
date: 11/2025

*/

void branch_and_cut(schedule_struct schedule, DataFrame simfile) {}

void branch_and_bound() {}
