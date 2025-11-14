#include "Nbody.hpp"
#include <algorithm>
#include <armadillo>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

double calculate_edelbaum_deltaV(double v0, double v1,
                                 double plane_diff_angle) {

  // assumes the transfer is between two circular orbits
  // assumes that there is no RAAN difference

  double deltaV =
      std::sqrt(pow(v0, 2) + pow(v1, 2) -
                2 * v1 * v0 * std::cos(M_PI / 2 * plane_diff_angle));

  return deltaV;
}

void low_thrust_di_carlo_independent() {};

void low_thrust_di_carlo_twophase() {};
