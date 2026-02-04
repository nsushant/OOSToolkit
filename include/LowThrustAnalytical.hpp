
#pragma once

#include <armadillo>
#include <string>
#include "Nbody.hpp"

double calculate_plane_diff_angle(double i_0, double i_1, double RAAN_0,double RAAN_1);

double calculate_edelbaum_deltaV(arma::vec v0, arma::vec v1, arma::vec r0, arma::vec r1, double t = 0, std::string method = "all");

// New orbital element-based Edelbaum function for enhanced performance
double calculate_edelbaum_deltaV_orbital_elems(
    orbital_elements orb_departure,  // departure satellite AT ARRIVAL TIME
    orbital_elements orb_arrival,    // arrival satellite AT ARRIVAL TIME  
    double transfer_time,           // t_arrival - t_departure
    std::string method
);

void low_thrust_di_carlo_independent();


void low_thrust_di_carlo_twophase();


