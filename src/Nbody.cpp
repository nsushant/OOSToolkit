
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include "Nbody.hpp"
#include "LambertSolver.hpp"

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


                /*
                std::cout << std::fixed << std::setprecision(3)
                     << t << " " << i << " " << sats[i].satname << " "
                     << sats[i].r(0) << " " << sats[i].r(1) << " " << sats[i].r(2) << " "
                     << sats[i].v(0) << " " << sats[i].v(1) << " " << sats[i].v(2) << "\n";
                */
            
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


    // testing lambert solver
    double tof = 120;

    // service satellite at t=0 
    arma::vec r1 = {5.63576e+06,4.2052e+06,915068};

    // client satellite at t=120
    arma::vec r2 = {7.02093e+06,540482,717244};

    std::vector<arma::vec> lambert_solutions = lambert_solver(r1, r2, tof, 3.986004418e14 , 0, 1);

    std::cout << "Number of solutions found: " << lambert_solutions.size() << std::endl ;

    for (int sols=0 ; sols < lambert_solutions.size(); sols++){
    
    std::cout << "Solution number: " << sols+1 << std::endl;
    std::cout << std::setprecision(10) << "V1 (at service station) (m/s): " << lambert_solutions[sols].subvec(0,2) << std::endl;
    std::cout << std::setprecision(10) << "V2 (at client) (m/s): " << lambert_solutions[sols].subvec(3,5) << std::endl;

    }



    return 0;
}
