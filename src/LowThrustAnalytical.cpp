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


double calculate_plane_diff_angle(double i_0, double i_1, double RAAN_0,double RAAN_1){

    double delta_i = i_1 - i_0; 
    double i_bar = (i_1 + i_0)/2.0; 
    double raan_change = RAAN_1 - RAAN_0; 
      
    //std::cout<<"RAAN Diff: "<<raan_change<<"\n";
    //std::cout<<"I diff: "<< delta_i << "\n";

    double g = std::sqrt(std::pow(delta_i,2) + std::sin(i_bar) * std::sin(i_bar) * std::pow(raan_change,2));


    return g; 


}


double calculate_edelbaum_deltaV(arma::vec v0, arma::vec v1, arma::vec r0, arma::vec r1) {

  // assumes the transfer is between two circular orbits
  // assumes that there is no RAAN difference
  

    orbital_elements orb1 = orb_elems_from_rv(r1, v1); 

    orbital_elements orb2 = orb_elems_from_rv(r0, v0); 

    double plane_diff_angle = calculate_plane_diff_angle(orb1.inclination, orb2.inclination, orb1.RAAN, orb2.RAAN);
    
    double V_0 = std::sqrt(MU_EARTH/orb1.semi_major_axis);
    double V_1 = std::sqrt(MU_EARTH/orb2.semi_major_axis); 


  double deltaV = std::sqrt(pow(V_0, 2) + pow(V_1, 2) - 2 * V_0 * V_1 * std::cos(M_PI / 2 * plane_diff_angle));

  //std::cout << "plane diff angle :" << plane_diff_angle<<"  ";

  return deltaV;
}

void low_thrust_di_carlo_independent() {};

void low_thrust_di_carlo_twophase() {};
