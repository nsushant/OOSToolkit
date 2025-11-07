#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include <chrono>
#include "Nbody.hpp"

/*
Author : S. Nigudkar (2025)

Functions that allow the initialization and propagation of a satellite constellation and a
service station in circular orbit. Available constellation arrangements are currently limited 
to - Delta Walker. 

Units used - m,s,kg

Time Integrator used - Runge-Kutta 4

*/

// sim functions // 

double deg_to_rads(double deg){

    //converts measured angles from degrees to radians

    return deg * M_PI / 180.0; 

}


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
                                             double sat_mass) {                       
                                                
    std::vector<satellite_object> sats;

    // predecided size of array
    sats.reserve(total_satnum);
    

    double a = altitude_m;
    // num satellites per plane
    int S = total_satnum / num_planes ;
    // the right accession (a total of 360 deg or 2pi rads) is divided up into n orbital planes. 
    double RAAN_step = 2.0 * M_PI / num_planes;
    
    // for each plane, 
    for (int p = 0; p < num_planes; ++p) {

        double RAAN = p * RAAN_step;

        // for each sat in this plane
        for (int s = 0; s < S; ++s) {

            // apply phasing
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


double calculate_edelbaum_deltaV(double v0,double v1, double plane_diff_angle){

    double deltaV = pow(v0,2) + pow(v1,2) - 2 * v1 * v0 * (M_PI/2 * plane_diff_angle);

    return sqrt(deltaV);

}



std::vector<satellite_object> circular_orbits(int total_satnum, 
                                             double altitude_m,
                                             double inclination_rad,
                                             double RAAN , double satmass){


    std::vector<satellite_object> sats;
    // predecided size of array
    sats.reserve(total_satnum);
                     
    double a = altitude_m;

    orbital_elements el;

    
    el.semi_major_axis = a;
    el.eccentricity = 0.0;
    el.inclination = inclination_rad;
    el.RAAN = RAAN;
    
    double phase = 2 * M_PI / total_satnum ; 

    el.augment_of_periapsis= 0.0;
    
    for (int p = 0; p < total_satnum; p++){

        el.true_anomaly = phase * p;

        std::string name = "sat_"+std::to_string(p) ; 
    
        sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords(name, el, satmass));
    
    }
                                                

    return sats;

}




void showProgressBar(int progress, int total, int barWidth ) {
    float ratio = static_cast<float>(progress) / total;
    int pos = static_cast<int>(barWidth * ratio);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "█";   // filled part
        else std::cout << "░";           // empty part
    }
    std::cout << "] " << std::setw(3) << int(ratio * 100) << "%\r";
    std::cout.flush();
    std::cout<<"\n";
}







void run_simulation(    std::string save_to_file, std::string arrangement, double t_final,double dt, 
                        double altitude_m,double num_planes, double num_satellites, 
                        double relative_phase, double inclination_in_radians, double satmass){


    double altitude = R_EARTH + altitude_m;
    double period = 2*M_PI * sqrt(pow(altitude,3) / MU_EARTH);

    std::cout << "orbital period: " <<  period;  

    force_model force_options;
    force_options.includeJ2 = false;
    force_options.includeMutual = false;

    std::vector<satellite_object> sats; 

    if (arrangement == "walker_delta"){

        sats = build_walker_constellation(num_planes, num_satellites, relative_phase,
                                                                        altitude,
                                                                        56.0 * M_PI/180.0,
                                                                        satmass);
        
    } 

    if (arrangement == "circular_orbit"){


        std::vector<satellite_object> sats = circular_orbits(num_satellites, 
                                             altitude,inclination_in_radians,
                                             0.0, satmass);

                                            
   
    }

    if (arrangement == "flower_constellation"){


        
        
        //std::vector<satellite_object> sats = build_coplanar_circular_orbit(); 

    }

    if (arrangement == "spiral_constellation"){
                 
        //std::vector<satellite_object> sats = build_coplanar_circular_orbit(); 

        
    }

    if (arrangement == "custom"){
        
        //std::vector<satellite_object> sats = build_coplanar_circular_orbit(); 

    }


    //init service satellite in circular orbit

    orbital_elements svc;
    svc.semi_major_axis = R_EARTH + 720000.0;
    svc.eccentricity = 0.001;
    svc.inclination = 48.0 * M_PI/180.0;
    svc.RAAN = 30.0 * M_PI/180.0;
    svc.augment_of_periapsis = 10.0 * M_PI/180.0;
    svc.true_anomaly = 0.0;
    sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords("service_1", svc, 500.0));

    //std::cout << "# time(s) index name x y z vx vy vz\n";

    //write header
    std::ofstream csv("../data/"+save_to_file);
    csv << "time_s,index,name,x,y,z,vx,vy,vz\n";


    double t = 0.0;
    int step = 0;


    std::cout<<"Running Simulation"; 
    while(t < (t_final - 1e-9)){
        // only write outputs to csv files every 10 timesteps
        // I am just thinning the datafile here so we can store it easily 

        if(step % 10 == 0){

            for(size_t i=0;i<sats.size();i++){

            
                csv << t << "," << i << "," << sats[i].satname << ","
                    << sats[i].r(0) << "," << sats[i].r(1) << "," << sats[i].r(2) << ","
                    << sats[i].v(0) << "," << sats[i].v(1) << "," << sats[i].v(2) << "\n";
            
            }

        }

        runge_kutta_step(sats, dt, force_options);
        
        t += dt;
        step++;
        //showProgressBar(t, t_final);
        

    }
    
    std::cout << std::endl << "Done!" << std::endl;

    csv.close();
    std::cout << "CSV saved to data/WalkerDelta.csv\n";



}








