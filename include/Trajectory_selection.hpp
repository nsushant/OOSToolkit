
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <armadillo>
#include "LambertSolver.hpp"
#include "data_access_lib.hpp"
#include "Nbody.hpp"
#include "Local_search.hpp"

/*
Author : S. Nigudkar (2025)


The aim of the below function is to choose lambert parameters r1,r2,tof and compute the associated delta Vs
for the most efficient orbital transfers given a requested service deadline t_service, the start time
t_prev, the name of the client satellite (departure destination) and that of the service satellite (arrival destination).

*/

void find_optimal_trajectory(std::string service_satname, std::string client_satname, double t_prev, double t_request, arma::vec &v1sol, arma::vec &v2sol, arma::vec &r1sol, arma::vec &r2sol, double &tof_optimal, std::vector<arma::vec> &trajs, DataFrame simfile, std::string method, bool write_to_file, double &DeltaVMinima);

double compact_optimal_calc(std::string satname1, std::string satname2, double departure_time, double arrival_time, DataFrame simfile);

void find_optimal_trajectory_no_iter(std::string service_satname, std::string client_satname, double t_departure, double t_arrival, DataFrame simfile, double &DeltaVMinima);

void run_exhaustive_search(std::string sat_from, std::string sat_to, double t_from, double t_to, double &deltaV_change, std::string simfilename, std::string outputfilename = "", std::string method = "lambert");
