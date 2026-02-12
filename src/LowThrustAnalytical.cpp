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

double raan_rate(double a, double e, double i) {

    return -3.0/2.0 * std::sqrt(MU_EARTH/std::pow(a,3)) * ( J2 * std::cos(i)) / std::pow((1.0 - std::pow(e,2)),2) * std::pow((R_EARTH/a),2) ;

};



double calculate_edelbaum_deltaV(arma::vec v0, arma::vec v1, arma::vec r0, arma::vec r1, double t,std::string method) {

  //the solutions provided cannot work beyond the following bounds


  // delta i > 10deg or
  // delta RAAN > 20-30 deg
  // delta a/a  > 0.2
  // k values > 3-5



  // calculate orbital elements

  orbital_elements orb1 = orb_elems_from_rv(r0, v0);

  orbital_elements orb2 = orb_elems_from_rv(r1, v1);

  double deltaV_tot = 0.0;


  double diff_raan_rate = (raan_rate(orb1.semi_major_axis,orb1.eccentricity,orb1.inclination) -
                          raan_rate(orb2.semi_major_axis,orb2.eccentricity,orb2.inclination));




  // is it feasible to make up the raan diff passively
  if (method == "passive"){


    double feasibility_with_nat_precession =  std::abs(orb1.RAAN - orb2.RAAN)/std::abs(diff_raan_rate);


    if (feasibility_with_nat_precession > t){

        method = "active";
    }

  }


  if (method == "No J2"){


    double plane_diff_angle = calculate_plane_diff_angle(orb1.inclination, orb2.inclination, orb1.RAAN, orb2.RAAN);

    double V_0 = std::sqrt(MU_EARTH/orb1.semi_major_axis);
    double V_1 = std::sqrt(MU_EARTH/orb2.semi_major_axis);


    double deltaV = std::sqrt(pow(V_0, 2) + pow(V_1, 2) - 2 * V_0 * V_1 * std::cos(M_PI / 2 * plane_diff_angle));


    deltaV_tot = deltaV;

  }


  if (method == "Debris Trajectory"){


    // the paper titled " Simple ΔV Approximation for Optimization of Debris-to-Debris Transfers "
    // by Shen and Casalino (2020) : arXiv:2004.02225 , proposes an approximate analytical law for
    // the quick calculation of delta V in the case of debris to debris transfer (e,i,a, RAAN all change quickly)


    //std::cout<< "i " << ", " << "RAAN" <<" , " << " e " << "\n";

    //std::cout<< orb1.inclination << ", " << orb1.RAAN << ", " << orb1.eccentricity << "\n";

    //std::cout<< orb2.inclination << ", " << orb2.RAAN << ", " << orb2.eccentricity << "\n";

    double i_0 = (orb2.inclination + orb1.inclination)/2.0;

    double a_0 = (orb2.semi_major_axis + orb1.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);



    double RAAN_dot_0 = (raan_rate(orb1.semi_major_axis,orb1.eccentricity,orb1.inclination) +
                        raan_rate(orb2.semi_major_axis,orb2.eccentricity,orb2.inclination)) / 2.0 ;



    double m = (7 * RAAN_dot_0)* std::sin(i_0) * t;

    double n = (RAAN_dot_0 * std::tan(i_0)) * t;


    // these orbital elements present a coordinate (non-physical) singularity
    // thus, when the required change is 0, all dependent quantities must be
    // manually set to 0

    //required RAAN change
    //In this form we actively perform the RAAN change using thrust
    //rather than taking advantage of the natural precession of the plane
    //due to J2 (gravitational torque due to the bulge of the Earth)

    double x  = (orb2.RAAN - orb1.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb2.semi_major_axis - orb1.semi_major_axis)/(2*a_0) * v_0;

    //required inclination change
    double z = (orb2.inclination - orb1.inclination)*v_0;

    // we consider three deltaV contributions
    // that deal with changing a,i,e,RAAN in a single transfer

    double deltaVa , deltaVb, deltaV_tot;

    double s_x, s_y, s_z;

    if ( x*x > 0.0 ){

      s_x = ( 2.0*x + m*y + n*z )/((4 + m*m + n*n)*x);

    }

    else{
      s_x = 0.0;
    }
    //std::cout << m << " " << n << "\n";

    if (y*y > 0.0){

      s_y = (2.0 * m * x - (4 + n*n) * y + m*n*z) / ((8 + 2*m*m + 2*n*n)*y);
    }

    else{

      s_y = 0.0;
    }

    if(z*z > 0.0){

       s_z = (2*n*x + (m*n*y) - (4 + m*m)*z)/((8 + (2*m*m) + (2*n*n))*z);

    }

    else{

        s_z = 0.0;
    }

    //std::cout << s_x <<" "<< s_y <<" "<<s_z << " "<< x <<" "<< y <<" " << z <<" " << "\n" ;


    double delta_x = m*s_y*y + n*s_z*z;

    deltaVa = std::sqrt(((s_x * x)*(s_x * x)) + ((s_y * y)*(s_y * y)) + ((s_z*z)*(s_z*z)));

    double deltaVbsq = ((x- (s_x * x) - delta_x)* (x- (s_x * x) - delta_x))  + ((y + s_y * y)* (y + s_y * y))+ ((z + s_z * z) *  (z + s_z * z)) ;

    deltaVb = std::sqrt(deltaVbsq);


    //small eccentricity changes

    double delta_e_y = orb2.eccentricity*std::sin(orb2.augment_of_periapsis) -  orb1.eccentricity*std::sin(orb1.augment_of_periapsis);

    double delta_e_x = orb2.eccentricity*std::cos(orb2.augment_of_periapsis) -  orb1.eccentricity*std::cos(orb1.augment_of_periapsis);

    double deltaV_e = 0.5 * v_0 * std::sqrt((delta_e_y*delta_e_y) +(delta_e_x*delta_e_x));

    std::cout << "delta Ve :" << deltaV_e << "Delta va "<< deltaVa << "delta Vb " << deltaVb  << "delta Vb sq" << deltaVbsq <<"\n";

    deltaV_tot = std::sqrt((deltaVa*deltaVa)+((deltaV_e*0.5)*(deltaV_e*0.5))) + std::sqrt((deltaVb*deltaVb)+((deltaV_e*0.5)*(deltaV_e*0.5)));

    // Simple NaN check - return penalty if result is invalid
    if (std::isnan(deltaV_tot) || std::isinf(deltaV_tot)) {
        //std::cout << "Warning: Edelbaum deltaV calculation returned NaN/Inf, returning penalty value" << std::endl;
        return 1e10; // Large penalty value
    }

  }



  if(method == "passive"){


    // Assuming a thrust-coast-thrust transfer with RAAN differences being covered using
    // natural precession alone. The method would only work if the travel time allows for
    // enough time for the precession to be leveraged.


    // we can check that the provided time of flight is feasible by simply checking if the
    // constraint x + delta x is met.


    double RAAN_dot_0 = (raan_rate(orb1.semi_major_axis,orb1.eccentricity,orb1.inclination) +
                        raan_rate(orb2.semi_major_axis,orb2.eccentricity,orb2.inclination)) / 2.0 ;



    // actively corrects the RAAN difference using thrusters

    double i_0 = (orb2.inclination + orb1.inclination)/2.0;
    double a_0 = (orb2.semi_major_axis + orb1.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);


    double x  = M_PI/2.0 * (orb2.RAAN - orb1.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb2.semi_major_axis - orb1.semi_major_axis)/(2.0*a_0) * v_0;

    //required inclination change
    double z = M_PI/2.0*(orb2.inclination - orb1.inclination)*v_0;


    // if the average inclination is 0 m=n=0 ,
    // the RAAN precession rate or omega_dot can only ever be 0 when i = 90
    // (so cos (i) is 0)

    double m,n;

    if ( (std::abs(i_0) < 1e-30) || (std::abs(RAAN_dot_0) < 1e-30) ){


      m = 0.0;
      n = 0.0;

    }

    else{

      m = M_PI/2.0 * (7 * RAAN_dot_0)* std::sin(i_0) * t;

      n = (RAAN_dot_0 * std::tan(i_0)) * std::sin(i_0) * t;

    }

    double m_n_denom = (m*m + n*n);

    double ky, kz;

    if ( std::abs(m_n_denom) > 0.0 ){


      // percent change in semi-major axis (first phase)
      ky = (y*y > 0.0) ? -(2*m*x - n*n*y + m*n*z) / (2*(m*m + n*n)*y) : 0.0;

      // percent change in inclination (first phase)
      kz = (z*z > 0.0) ? -(2*n*x + m*n*y - m*m*z) / (2*(m*m + n*n )*z) : 0.0;

    }

    else{

      ky = 0.0;
      kz = 0.0;
    }

    double delta_x = m*ky*y + n*kz*z;


    double raan_check = x + delta_x;

    //if(raan_check > 1e-12){

      //std::cout<<" raan not matched - check on line 247 " << "\n";
      //std::cout<<" not enough time for RAAN to be matched by natural precession "<<"\n";

    //}

    // semi major axis and inclination change (first phase)
    double deltav_1 = std::sqrt((ky*y)*(ky*y) + (kz*z)*(kz*z) );


    // remaining change to match target orbit (second phase)
    double deltav_2 = std::sqrt((y-ky*y)*(y-ky*y) + (z-kz*z)*(z-kz*z));


    deltaV_tot = deltav_1 + deltav_2;


  }

  if(method == "active"){

    // here we assume that the RAAN is matched by using the thrusters
    // and natural precession does not provide enough of a RAAN change (shorter tof)


    double RAAN_dot_0 = (raan_rate(orb1.semi_major_axis,orb1.eccentricity,orb1.inclination) +
                        raan_rate(orb2.semi_major_axis,orb2.eccentricity,orb2.inclination)) / 2.0 ;


    double i_0 = (orb2.inclination + orb1.inclination)/2.0;
    double a_0 = (orb2.semi_major_axis + orb1.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);

    double kx_x, ky_y, kz_z;


    // required active raan change
    double x  = M_PI/2.0 * (orb2.RAAN - orb1.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb2.semi_major_axis - orb1.semi_major_axis)/(2.0*a_0) * v_0;

    //required inclination change
    double z = M_PI/2.0*(orb2.inclination - orb1.inclination)*v_0;


    double m,n;



    if ( (std::abs(i_0) < 1e-30) || (std::abs(RAAN_dot_0) < 1e-30) )
    {

        m = 0.0;
        n= 0.0;

    }

    else
    {
        m = M_PI/2.0 * (7 * RAAN_dot_0)* std::sin(i_0) * t;

        n = (RAAN_dot_0 * std::tan(i_0)) * std::sin(i_0) * t;
    }


    kx_x = (2.0*x + m*y + n*z)/(4 + m*m + n*n);
    ky_y = - (2*m*x - (4 + n*n)*y + m*n*z)/(8+2*m*m + 2*n*n);
    kz_z = - (2*n*x + m*n*y - (4 + m*m)*z)/(8 + 2*m*m + 2*n*n);



    double delta_x = m*ky_y + n*kz_z;


    double deltav1 = std::sqrt((kx_x*kx_x) + (ky_y*ky_y) + (kz_z*kz_z));

    double deltav2 = std::sqrt((x - kx_x + delta_x)*(x - kx_x + delta_x) + (y-ky_y)*(y-ky_y) + (z-kz_z)*(z-kz_z));

    deltaV_tot = deltav1 + deltav2;


  }




  return deltaV_tot;

}




void low_thrust_di_carlo_independent() {};




void low_thrust_di_carlo_twophase() {};

// New orbital element-based Edelbaum function for enhanced performance
double calculate_edelbaum_deltaV_orbital_elems(
    orbital_elements orb_departure,  // departure satellite AT ARRIVAL TIME
    orbital_elements orb_arrival,    // arrival satellite AT ARRIVAL TIME
    double transfer_time,           // t_arrival - t_departure
    std::string method
) {
  //the solutions provided cannot work beyond the following bounds

  // delta i > 10deg or
  // delta RAAN > 20-30 deg
  // delta a/a  > 0.2
  // k values > 3-5

  double deltaV_tot = 0.0;
  auto raan_rate = [] (double a, double e, double i) {
      return -3.0/2.0 * std::sqrt(MU_EARTH/std::pow(a,3)) * ( J2 * std::cos(i)) / std::pow((1.0 - std::pow(e,2)),2) * std::pow((R_EARTH/a),2) ;
  };

  double diff_raan_rate = (raan_rate(orb_departure.semi_major_axis,orb_departure.eccentricity,orb_departure.inclination) -
                          raan_rate(orb_arrival.semi_major_axis,orb_arrival.eccentricity,orb_arrival.inclination));



  // is it feasible to make up the raan diff passively
  if (method == "passive"){

    double feasibility_with_nat_precession =  std::abs(orb_departure.RAAN - orb_arrival.RAAN)/std::abs(diff_raan_rate);

    if (feasibility_with_nat_precession > transfer_time){

        method = "active";
    }

  }

  if (method == "No J2"){

    double plane_diff_angle = calculate_plane_diff_angle(orb_departure.inclination, orb_arrival.inclination, orb_departure.RAAN, orb_arrival.RAAN);

    double V_0 = std::sqrt(MU_EARTH/orb_departure.semi_major_axis);
    double V_1 = std::sqrt(MU_EARTH/orb_arrival.semi_major_axis);


    double deltaV = std::sqrt(pow(V_0, 2) + pow(V_1, 2) - 2 * V_0 * V_1 * std::cos(M_PI / 2 * plane_diff_angle));


    deltaV_tot = deltaV;

  }


  if (method == "Debris Trajectory"){

    // the paper titled " Simple ΔV Approximation for Optimization of Debris-to-Debris Transfers "
    // by Shen and Casalino (2020) : arXiv:2004.02225 , proposes an approximate analytical law for
    // the quick calculation of delta V in the case of debris to debris transfer (e,i,a, RAAN all change quickly)

    double i_0 = (orb_arrival.inclination + orb_departure.inclination)/2.0;

    double a_0 = (orb_arrival.semi_major_axis + orb_departure.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);



    double RAAN_dot_0 = (raan_rate(orb_departure.semi_major_axis,orb_departure.eccentricity,orb_departure.inclination) +
                        raan_rate(orb_arrival.semi_major_axis,orb_arrival.eccentricity,orb_arrival.inclination)) / 2.0 ;



    double m = (7 * RAAN_dot_0)* std::sin(i_0) * transfer_time;

    double n = (RAAN_dot_0 * std::tan(i_0)) * transfer_time;


    // these orbital elements present a coordinate (non-physical) singularity
    // thus, when the required change is 0, all dependent quantities must be
    // manually set to 0

    //required RAAN change
    //In this form we actively perform the RAAN change using thrust
    //rather than taking advantage of the natural precession of the plane
    //due to J2 (gravitational torque due to the bulge of the Earth)

    double x  = (orb_arrival.RAAN - orb_departure.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb_arrival.semi_major_axis - orb_departure.semi_major_axis)/(2*a_0) * v_0;

    //required inclination change
    double z = (orb_arrival.inclination - orb_departure.inclination)*v_0;

    // we consider three deltaV contributions
    // that deal with changing a,i,e,RAAN in a single transfer

    double deltaVa , deltaVb, deltaV_tot;

    double s_x, s_y, s_z;

    if ( x*x > 0.0 ){

      s_x = ( 2.0*x + m*y + n*z )/((4 + m*m + n*n)*x);

    }

    else{
      s_x = 0.0;
    }
    //std::cout << m << " " << n << "\n";

    if (y*y > 0.0){

      s_y = (2.0 * m * x - (4 + n*n) * y + m*n*z) / ((8 + 2*m*m + 2*n*n)*y);
    }

    else{

      s_y = 0.0;
    }

    if(z*z > 0.0){

       s_z = (2*n*x + (m*n*y) - (4 + m*m)*z)/((8 + (2*m*m) + (2*n*n))*z);

    }

    else{

        s_z = 0.0;
    }

    //std::cout << s_x <<" "<< s_y <<" "<<s_z << " "<< x <<" "<< y <<" " << z <<" " << "\n" ;


    double delta_x = m*s_y*y + n*s_z*z;

    deltaVa = std::sqrt(((s_x * x)*(s_x * x)) + ((s_y * y)*(s_y * y)) + ((s_z*z)*(s_z*z)));

    double deltaVbsq = ((x- (s_x * x) - delta_x)* (x- (s_x * x) - delta_x))  + ((y + s_y * y)* (y + s_y * y))+ ((z + s_z * z) *  (z + s_z * z)) ;

    deltaVb = std::sqrt(deltaVbsq);


    //small eccentricity changes

    double delta_e_y = orb_arrival.eccentricity*std::sin(orb_arrival.augment_of_periapsis) -  orb_departure.eccentricity*std::sin(orb_departure.augment_of_periapsis);

    double delta_e_x = orb_arrival.eccentricity*std::cos(orb_arrival.augment_of_periapsis) -  orb_departure.eccentricity*std::cos(orb_departure.augment_of_periapsis);

    double deltaV_e = 0.5 * v_0 * std::sqrt((delta_e_y*delta_e_y) +(delta_e_x*delta_e_x));

    std::cout << "delta Ve :" << deltaV_e << "Delta va "<< deltaVa << "delta Vb " << deltaVb  << "delta Vb sq" << deltaVbsq <<"\n";

    deltaV_tot = std::sqrt((deltaVa*deltaVa)+((deltaV_e*0.5)*(deltaV_e*0.5))) + std::sqrt((deltaVb*deltaVb)+((deltaV_e*0.5)*(deltaV_e*0.5)));

    // Simple NaN check - return penalty if result is invalid
    if (std::isnan(deltaV_tot) || std::isinf(deltaV_tot)) {
        //std::cout << "Warning: Edelbaum deltaV calculation returned NaN/Inf, returning penalty value" << std::endl;
        return 1e10; // Large penalty value
    }

  }



  if(method == "passive"){


    // Assuming a thrust-coast-thrust transfer with RAAN differences being covered using
    // natural precession alone. The method would only work if the travel time allows for
    // enough time for the precession to be leveraged.


    // we can check that the provided time of flight is feasible by simply checking if the
    // constraint x + delta x is met.


    double RAAN_dot_0 = (raan_rate(orb_departure.semi_major_axis,orb_departure.eccentricity,orb_departure.inclination) +
                        raan_rate(orb_arrival.semi_major_axis,orb_arrival.eccentricity,orb_arrival.inclination)) / 2.0 ;



    // actively corrects the RAAN difference using thrusters

    double i_0 = (orb_arrival.inclination + orb_departure.inclination)/2.0;
    double a_0 = (orb_arrival.semi_major_axis + orb_departure.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);


    double x  = M_PI/2.0 * (orb_arrival.RAAN - orb_departure.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb_arrival.semi_major_axis - orb_departure.semi_major_axis)/(2.0*a_0) * v_0;

    //required inclination change
    double z = M_PI/2.0*(orb_arrival.inclination - orb_departure.inclination)*v_0;


    // if the average inclination is 0 m=n=0 ,
    // the RAAN precession rate or omega_dot can only ever be 0 when i = 90
    // (so cos (i) is 0)

    double m,n;

    if ( (std::abs(i_0) < 1e-30) || (std::abs(RAAN_dot_0) < 1e-30) ){


      m = 0.0;
      n = 0.0;

    }

    else{

      m = M_PI/2.0 * (7 * RAAN_dot_0)* std::sin(i_0) * transfer_time;

      n = (RAAN_dot_0 * std::tan(i_0)) * std::sin(i_0) * transfer_time;

    }

    double m_n_denom = (m*m + n*n);

    double ky, kz;

    if ( std::abs(m_n_denom) > 0.0 ){


      // percent change in semi-major axis (first phase)
      ky = (y*y > 0.0) ? -(2*m*x - n*n*y + m*n*z) / (2*(m*m + n*n)*y) : 0.0;

      // percent change in inclination (first phase)
      kz = (z*z > 0.0) ? -(2*n*x + m*n*y - m*m*z) / (2*(m*m + n*n )*z) : 0.0;

    }

    else{

      ky = 0.0;
      kz = 0.0;
    }

    double delta_x = m*ky*y + n*kz*z;


    double raan_check = x + delta_x;

    //if(raan_check > 1e-12){

      //std::cout<<" raan not matched - check on line 247 " << "\n";
      //std::cout<<" not enough time for RAAN to be matched by natural precession "<<"\n";

    //}

    // semi major axis and inclination change (first phase)
    double deltav_1 = std::sqrt((ky*y)*(ky*y) + (kz*z)*(kz*z) );


    // remaining change to match target orbit (second phase)
    double deltav_2 = std::sqrt((y-ky*y)*(y-ky*y) + (z-kz*z)*(z-kz*z));


    deltaV_tot = deltav_1 + deltav_2;

  }

  if(method == "active"){

    // here we assume that the RAAN is matched by using the thrusters
    // and natural precession does not provide enough of a RAAN change (shorter tof)


    double RAAN_dot_0 = (raan_rate(orb_departure.semi_major_axis,orb_departure.eccentricity,orb_departure.inclination) +
                        raan_rate(orb_arrival.semi_major_axis,orb_arrival.eccentricity,orb_arrival.inclination)) / 2.0 ;


    double i_0 = (orb_arrival.inclination + orb_departure.inclination)/2.0;
    double a_0 = (orb_arrival.semi_major_axis + orb_departure.semi_major_axis)/2.0;

    double v_0 = std::sqrt(MU_EARTH/a_0);

    double kx_x, ky_y, kz_z;



    // required active raan change
    double x  = M_PI/2.0 * (orb_arrival.RAAN - orb_departure.RAAN) * std::sin(i_0) * v_0;

    //required semi-maj axis change
    double y =  (orb_arrival.semi_major_axis - orb_departure.semi_major_axis)/(2.0*a_0) * v_0;

    //required inclination change
    double z = M_PI/2.0*(orb_arrival.inclination - orb_departure.inclination)*v_0;


    double m,n;



    if ( (std::abs(i_0) < 1e-30) || (std::abs(RAAN_dot_0) < 1e-30) )
    {

        m = 0.0;
        n= 0.0;

    }

    else
    {
        m = M_PI/2.0 * (7 * RAAN_dot_0)* std::sin(i_0) * transfer_time;

        n = (RAAN_dot_0 * std::tan(i_0)) * std::sin(i_0) * transfer_time;
    }



    kx_x = (2.0*x + m*y + n*z)/(4 + m*m + n*n);
    ky_y = - (2*m*x - (4 + n*n)*y + m*n*z)/(8+2*m*m + 2*n*n);
    kz_z = - (2*n*x + m*n*y - (4 + m*m)*z)/(8 + 2*m*m + 2*n*n);



    double delta_x = m*ky_y + n*kz_z;


    double deltav1 = std::sqrt((kx_x*kx_x) + (ky_y*ky_y) + (kz_z*kz_z));

    double deltav2 = std::sqrt((x - kx_x + delta_x)*(x - kx_x + delta_x) + (y-ky_y)*(y-ky_y) + (z-kz_z)*(z-kz_z));

    deltaV_tot = deltav1 + deltav2;

  }

  return deltaV_tot;
}
