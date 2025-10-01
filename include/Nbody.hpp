

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>

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

arma::vec3 acceleration_due_to_central_body (const arma::vec3& r){
    
    double r_magnitude = arma::norm(r);
    
    if(r_magnitude == 0){ 

        return arma::vec3();
    }

    return r * (-MU_EARTH / (pow(r_magnitude,3)));
}

// acceleration taking j2 perturbations into account
arma::vec3 acceleration_due_to_J2(const arma::vec3& r){

    double rx = r(0);
    double ry = r(1);
    double rz = r(2);

    double r2 = rx*rx + ry*ry + rz*rz;
    double r5 = pow(r2, 2.5);

    double z2 = rz*rz;
    
    double factor = 1.5 * J2 * MU_EARTH * R_EARTH*R_EARTH / r5;
    
    double common = 5.0 * z2 / r2 - 1.0;
    
    arma::vec3 a;
    
    a(0) = rx * factor * common;
    a(1)= ry * factor * common;
    a(2) = rz * factor * (5.0 * z2 / r2 - 3.0);
    return a;
}


// compute accelerations given froce options
arma::vec3 compute_acceleration(const satellite_object& sat, const std::vector<satellite_object>& all_sats, const force_model& force_options){
    
    arma::vec3 acceleration_total = acceleration_due_to_central_body(sat.r);

    // adding a J2 perturbation
    // (Non-keplerian)

    if(force_options.includeJ2){
        acceleration_total += acceleration_due_to_J2(sat.r);
    }


    // the perturbations in orbits caused by the mutual attraction between satellites
    // (Non keplerian)
    if(force_options.includeMutual){

        // calculate the mutual force applied by all sats
        for(const auto& other_satellite : all_sats){
            // exclude self from calculation (force applied by satellite on itself)
            if(&other_satellite == &sat) continue;

            // distance between satellites
            arma::vec3 dr_sat = other_satellite.r - sat.r;

            double d = arma::norm(dr_sat);
            
            if(d == 0) continue;

                acceleration_total += dr_sat * (G_CONST * other_satellite.mass / std::pow(d,3));
        }
    }
    return acceleration_total;
}


// 4th order runge kutta integration to propogate satellite objects forwards in time. 
void runge_kutta_step(std::vector<satellite_object>& sats, double dt, const force_model& force_options){
    
    int num_sats = sats.size();

    // initialize vectors to store RK additions for both position and velocity
    std::vector<arma::vec> k1r(num_sats), k1v(num_sats), k2r(num_sats), k2v(num_sats), k3r(num_sats), k3v(num_sats), k4r(num_sats), k4v(num_sats);
    
    // we will store the outputs of the rk4 passes in this vector  
    std::vector<satellite_object> temp = sats;

    // iterating through all satellites 
    
    // RK4 Algo
    
    // first pass = f(yn)
    // k1 = f(yn + h * (fist pass) / 2 )
    // k2 = f(yn + h * (k1) / 2 )
    // k3 = f(yn + h * k2 )
    // k4 = y_(n+1) = y_n + h/6 (first pass + 2k1 + 2k2 + k3)

    // first pass is already stored in sats 

    // now we use the values in sats to calculate the rk1 
    for(int i=0;i<num_sats;i++) k1r[i] = sats[i].v;    
    for(int i=0;i<num_sats;i++) k1v[i] = compute_acceleration(sats[i], sats, force_options);
    for(int i=0;i<num_sats;i++){
        
        temp[i].r = sats[i].r + k1r[i]*(0.5*dt);
        temp[i].v = sats[i].v + k1v[i]*(0.5*dt);
    }

    // rk2 

    for(int i=0;i<num_sats;i++) k2r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k2v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){
        
        temp[i].r = sats[i].r + k2r[i]*(0.5*dt);
        temp[i].v = sats[i].v + k2v[i]*(0.5*dt);
    }


    // rk3 

    for(int i=0;i<num_sats;i++) k3r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k3v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){
        
        temp[i].r = sats[i].r + k3r[i]*dt;
        temp[i].v = sats[i].v + k3v[i]*dt;
    }

    
    // using above passes to compute y_n+1

    for(int i=0;i<num_sats;i++) k4r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k4v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){

        sats[i].r = sats[i].r + (k1r[i] + k2r[i]*2.0 + k3r[i]*2.0 + k4r[i]) * (dt/6.0);
        sats[i].v = sats[i].v + (k1v[i] + k2v[i]*2.0 + k3v[i]*2.0 + k4v[i]) * (dt/6.0);
    }



}


// distributing satellites in a walker constellation
std::vector<satellite_object> build_walker_constellation(int num_planes, int total_satnum, int phase ,
                                             double altitude_m,
                                             double inclination_rad,
                                             double sat_mass=100.0) {                       
                                                
    std::vector<satellite_object> sats;

    // predecided size of array
    sats.reserve(total_satnum);
    

    double a = R_EARTH + altitude_m;
    int S = total_satnum / num_planes ;
    double RAAN_step = 2.0 * M_PI / num_planes;
    for (int p = 0; p < num_planes; ++p) {

        double RAAN = p * RAAN_step;

        for (int s = 0; s < S; ++s) {

            double phase = 2.0 * M_PI * ( (s + (double)phase * p / num_planes) / S );
            
            orbital_elements el;
            el.semi_major_axis = a;
            el.eccentricity = 0.0;
            el.inclination = inclination_rad;
            el.RAAN = RAAN;
            el.augment_of_periapsis= 0.0;
            el.true_anomaly = phase;
            std::string name = "sat_" + std::to_string(p*S + s);
            sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords(name, el, sat_mass));
        }
    }
    return sats;
}
