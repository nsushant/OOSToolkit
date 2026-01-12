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


double calculate_edelbaum_deltaV(arma::vec v0, arma::vec v1, arma::vec r0, arma::vec r1, double t,std::string method) {

  
  
  // assumes the transfer is between two circular orbits
  // assumes that there is no RAAN difference
  

  orbital_elements orb1 = orb_elems_from_rv(r1, v1); 

  orbital_elements orb2 = orb_elems_from_rv(r0, v0); 
  
  if (method == "No J2"){


    double plane_diff_angle = calculate_plane_diff_angle(orb1.inclination, orb2.inclination, orb1.RAAN, orb2.RAAN);
    
    double V_0 = std::sqrt(MU_EARTH/orb1.semi_major_axis);
    double V_1 = std::sqrt(MU_EARTH/orb2.semi_major_axis); 


  double deltaV = std::sqrt(pow(V_0, 2) + pow(V_1, 2) - 2 * V_0 * V_1 * std::cos(M_PI / 2 * plane_diff_angle));

  //std::cout << "plane diff angle :" << plane_diff_angle<<"  ";

  return deltaV;

  }


  else{



    double i_0 = (orb1.inclination + orb2.inclination)/2.0;
    
    double a_0 = (orb1.semi_major_axis + orb2.semi_major_axis)/2.0; 

    double v_0 = std::sqrt(MU_EARTH/a_0); 

    auto raan_rate = [] (double a, double e, double i) {

      return -3.0/2.0 * std::sqrt(MU_EARTH/std::pow(a,3)) * ( J2 * std::cos(i)) / std::pow((1.0 - std::pow(e,2)),2) * std::pow((R_EARTH/a),2) ;



    };


    double RAAN_dot_0 = (raan_rate(orb1.semi_major_axis,orb1.eccentricity,orb1.inclination) + 
                        raan_rate(orb2.semi_major_axis,orb2.eccentricity,orb2.inclination)) / 2.0 ; 



    double m = (7 * RAAN_dot_0)* std::sin(i_0) * t; 

    double n = (RAAN_dot_0 * std::tan(i_0)) * t; 


    double x  = (orb2.RAAN - orb1.RAAN) * std::sin(i_0) * v_0;

    double y =  (orb2.semi_major_axis - orb1.semi_major_axis)/(2*a_0) * v_0;

    double z = (orb2.inclination - orb1.inclination)*v_0; 


    double deltaVa , deltaVb, deltaV_tot;



    double s_x = ( 2.0*x + m*y + n*z )/((4 + m*m + n*n)*x); 
    double s_y = (2.0 * m * x - (4 + n*n) * y + m*n*z) / (8 + 2*m*m + 2*n*n)*y; 
    double s_z = (2*n*x + m*n*y - (4 + m*m)*z)/((8 + 2*m*m + 2*n*n)*z);
  
    double delta_x = m*s_y*y + n*s_z*z;

    deltaVa = std::sqrt((s_x * x)*(s_x * x) + (s_y * y)*(s_y * y) + (s_z*z)*(s_z*z));

    deltaVb = std::sqrt( (x-s_x * x - delta_x)* (x-s_x * x - delta_x)  +  (y + s_y * y)*  (y + s_y * y) + (z + s_z * z) *  (z + s_z * z) );

    //small eccentricity changes 
    
    double delta_e_y = orb2.eccentricity*std::sin(orb2.augment_of_periapsis) -  orb1.eccentricity*std::sin(orb1.augment_of_periapsis); 

    double delta_e_x = orb2.eccentricity*std::cos(orb2.augment_of_periapsis) -  orb1.eccentricity*std::cos(orb1.augment_of_periapsis); 

    double deltaV_e = 0.5 * v_0 * std::sqrt(delta_e_y*delta_e_y +delta_e_x*delta_e_x);

    deltaV_tot = std::sqrt(deltaVa*deltaVa+(deltaV_e*0.5)*(deltaV_e*0.5)) + std::sqrt(deltaVb*deltaVb+(deltaV_e*0.5)*(deltaV_e*0.5));

  
    return deltaV_tot; 



  }


}




void low_thrust_di_carlo_independent() {};




void low_thrust_di_carlo_twophase() {};
