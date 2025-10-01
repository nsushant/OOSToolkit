
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include "Nbody.hpp"


int main(){
    double altitude_km = 700.0;
    int num_planes = 5;   
    int num_satellites = 50;  
    int relative_phase = 1;   

    double altitude = R_EARTH + altitude_km*1000.0;
    double period = 2*M_PI * sqrt(pow(altitude,3) / MU_EARTH);
    double dt = 10.0;
    double t_final = period * 1.05;

    force_model force_options;
    force_options.includeJ2 = false;
    force_options.includeMutual = false;

    std::vector<satellite_object> sats = build_walker_constellation(num_planes, num_satellites, relative_phase,
                                                        altitude_km * 1000.0,

                                                        53.0 * M_PI/180.0,

                                                        100.0);

    orbital_elements svc;
    svc.semi_major_axis = R_EARTH + 720000.0;
    svc.eccentricity = 0.001;
    svc.inclination = 48.0 * M_PI/180.0;
    svc.RAAN = 30.0 * M_PI/180.0;
    svc.augment_of_periapsis = 10.0 * M_PI/180.0;
    svc.true_anomaly = 0.0;
    sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords("service_1", svc, 500.0));

    std::cout << "# time(s) index name x y z vx vy vz\n";

    //write header
    std::ofstream csv("../data/WalkerDelta.csv");
    csv << "time_s,index,name,x,y,z,vx,vy,vz\n";


    double t = 0.0;
    int step = 0;


    while(t < (t_final - 1e-9)){

        // only write outputs to csv files every 6 timesteps
        // I am just thinning the datafile here so we can store it easily 

        if(step % 6 == 0){

            for(size_t i=0;i<sats.size();i++){

                std::cout << std::fixed << std::setprecision(3)
                     << t << " " << i << " " << sats[i].satname << " "
                     << sats[i].r(0) << " " << sats[i].r(1) << " " << sats[i].r(2) << " "
                     << sats[i].v(0) << " " << sats[i].v(1) << " " << sats[i].v(2) << "\n";
            
            
                csv << t << "," << i << "," << sats[i].satname << ","
                    << sats[i].r(0) << "," << sats[i].r(1) << "," << sats[i].r(2) << ","
                    << sats[i].v(0) << "," << sats[i].v(1) << "," << sats[i].v(2) << "\n";
            
            }

        }

        runge_kutta_step(sats, dt, force_options);
        
        t += dt;
        step++;


    }

    csv.close();
    std::cout << "CSV saved to data/WalkerDelta.csv\n";
    return 0;
}
