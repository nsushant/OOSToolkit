#pragma once 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include <chrono>


/*
Author : S. Nigudkar (2025)

Functions that allow the initialization and propagation of a satellite constellation and a
service station in circular orbit. Available constellation arrangements are currently limited 
to - Delta Walker. 

Units used - m,s,kg

Time Integrator used - Runge-Kutta 4

*/



// simulation constants // 
const double MU_EARTH = 3.986004418e14; // m^3/s^2
const double R_EARTH  = 6378137.0;      // m
const double J2       = 1.08262668e-3;
const double G_CONST  = 6.67430e-11;    // gravitational constant m^3/kg/s^2


// simulation structs //

struct orbital_elements
{   
    double inclination;
    double RAAN;
    double augment_of_periapsis;
    double true_anomaly;
    double semi_major_axis;
    double eccentricity;  
};

struct satellite_object
{
    std::string satname;
    double mass;      // kg
    arma::vec3 r;           // position in ECI (m)
    arma::vec3 v;           // velocity in ECI (m/s)

    //constructor
    satellite_object(std::string n="", double m=1.0): satname(n), mass(m), r(), v() {}

    // function to convert satellite coordinates from the perifocal coordinate system to the ECI frame
    // static makes the function available, independent from the instance of a satellite_object
    static satellite_object from_satellite_normal_to_ECI_coords(const std::string& satname, const orbital_elements& orb_elems, double mass=1.0){

        //mass and name remain the same 
        satellite_object sat_ECI(satname,mass); 
        double a = orb_elems.semi_major_axis;
        double e = orb_elems.eccentricity;
        double i = orb_elems.inclination;
        double RAAN = orb_elems.RAAN;
        double w = orb_elems.augment_of_periapsis;
        double nu = orb_elems.true_anomaly;

        // we can also calculate the semi-parameter p and the specific momentum h as well as the magnitude of the position vec
        // in the perifocal coordinate frame
     
        double p = a * (1 - e*e); 
        double mag_r = p / (1 + e * cos(nu)); // trajectory equation
        double h = sqrt(MU_EARTH * p);

        // statement of r and v in perifocal coordinate system 
        arma::vec r_perifocal = { mag_r * cos(nu), mag_r * sin(nu), 0.0 };
        arma::vec v_perifocal = { -MU_EARTH / h * sin(nu), MU_EARTH / h * (e + cos(nu)), 0.0 };

        // rotation matricies 
        arma::mat Rz_RAAN = { {cos(RAAN), -sin(RAAN), 0},
                              {sin(RAAN),  cos(RAAN), 0},
                              {0,          0,         1} };

        arma::mat Rx = { {1, 0, 0},
                         {0, cos(i), -sin(i)},
                         {0, sin(i),  cos(i)} };

        arma::mat Rz_w = {  {cos(w), -sin(w), 0},
                            {sin(w),  cos(w), 0},
                            {0,       0,      1} };

        // Rotation from perifocal to ECI
        arma::mat R_peri_to_ECI = Rz_RAAN * Rx * Rz_w;

        // Apply rotation 
        sat_ECI.r = R_peri_to_ECI * r_perifocal;
        sat_ECI.v = R_peri_to_ECI * v_perifocal;

        return sat_ECI;
        
    }

};

struct force_model {
    bool includeJ2 = true;
    bool includeMutual = true;
};


// sim functions // 

arma::vec3 acceleration_due_to_central_body (const arma::vec3& r);

// acceleration taking j2 perturbations into account
arma::vec3 acceleration_due_to_J2(const arma::vec3& r);


// compute accelerations given froce options
arma::vec3 compute_acceleration(const satellite_object& sat, 
                                const std::vector<satellite_object>& all_sats, 
                                const force_model& force_options    );

// 4th order runge kutta integration to propogate satellite objects forwards in time. 
void runge_kutta_step(std::vector<satellite_object>& sats, double dt, const force_model& force_options);

// distributing satellites in a walker constellation
std::vector<satellite_object> build_walker_constellation(   int num_planes, 
                                                            int total_satnum, 
                                                            int phase,
                                                            double altitude_m,
                                                            double inclination_rad,
                                                            double sat_mass=100.0   );                        
                                                


double calculate_edelbaum_deltaV(double v0,double v1, double plane_diff_angle);


void showProgressBar(int progress, int total, int barWidth = 50);

