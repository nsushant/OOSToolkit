
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include "LambertSolver.hpp"
#include "data_access_lib.hpp"




/*
Author : S. Nigudkar (2025)


The aim of the below function is to choose lambert parameters r1,r2,tof and compute the associated delta Vs 
for the most efficient orbital transfers given a requested service deadline t_service, the start time
t_prev, the name of the client satellite (departure destination) and that of the service satellite (arrival destination).

*/





void find_optimal_trajectory(std::string client_satname,std::string service_satname, double t_prev, double t_request ,arma::vec &v1sol, arma::vec &v2sol, arma::vec &r1sol, arma::vec &r2sol, double &tof_optimal, DataFrame simfile){


    // methods used on dataframes and the definition of the data structure
    // are both located in the data_access_lib.hpp file 
    std::vector<std::string> satnames = simfile["name"];


    // .getNumeric casts the strings stored in the csv file to floats
    arma::vec x_sat = simfile.getNumeric("x");
    arma::vec y_sat = simfile.getNumeric("y");
    arma::vec z_sat = simfile.getNumeric("z");
    arma::vec vx_sat = simfile.getNumeric("vx");
    arma::vec vy_sat = simfile.getNumeric("vy");
    arma::vec vz_sat = simfile.getNumeric("vz");
    arma::vec t_sat = simfile.getNumeric("time_s");


    // find_idxs_of_match finds the indicies at which satnames == service_satname
    arma::uvec service_sat_idxs = find_idxs_of_match(satnames,service_satname);
    arma::uvec client_sat_idxs = find_idxs_of_match(satnames,client_satname);

    // locations of service station
    arma::vec t_service_sat = t_sat.elem(service_sat_idxs);
    arma::uvec t_service_prev = arma::find(t_service_sat == t_prev);

    arma::vec x_service_sat = x_sat.elem(service_sat_idxs);
    arma::vec y_service_sat = y_sat.elem(service_sat_idxs);
    arma::vec z_service_sat = z_sat.elem(service_sat_idxs);


    // service station velocities
    arma::vec vx_service_sat = vx_sat.elem(service_sat_idxs);
    arma::vec vy_service_sat = vy_sat.elem(service_sat_idxs);
    arma::vec vz_service_sat = vz_sat.elem(service_sat_idxs);

    // service demand
    arma::vec t_client = t_sat.elem(client_sat_idxs) ; 
    arma::vec x_client = x_sat.elem(client_sat_idxs) ; 
    arma::vec y_client = y_sat.elem(client_sat_idxs) ; 
    arma::vec z_client = z_sat.elem(client_sat_idxs) ; 


    // client velocities
    arma::vec vx_client = vx_sat.elem(client_sat_idxs) ; 
    arma::vec vy_client = vy_sat.elem(client_sat_idxs) ; 
    arma::vec vz_client = vz_sat.elem(client_sat_idxs) ; 

    // available x,t for the service shuttle 
    arma::uvec t_service_idxs = arma::find( t_service_sat < t_request);
    arma::vec available_t_service = t_service_sat.elem(t_service_idxs);

    arma::vec available_x_service = x_service_sat.elem(t_service_idxs);
    arma::vec available_y_service = y_service_sat.elem(t_service_idxs);
    arma::vec available_z_service = z_service_sat.elem(t_service_idxs);

    arma::vec available_vx_service = vx_service_sat.elem(t_service_idxs);
    arma::vec available_vy_service = vy_service_sat.elem(t_service_idxs);
    arma::vec available_vz_service = vz_service_sat.elem(t_service_idxs);


    // available x,t for the client
    arma::uvec t_client_idxs = arma::find(t_client <= t_request);
    arma::vec available_t_client = t_client.elem(t_client_idxs);

    arma::vec available_x_client = x_client.elem(t_client_idxs);
    arma::vec available_y_client = y_client.elem(t_client_idxs);
    arma::vec available_z_client = z_client.elem(t_client_idxs);

    arma::vec available_vx_client = vx_client.elem(t_client_idxs);
    arma::vec available_vy_client = vy_client.elem(t_client_idxs);
    arma::vec available_vz_client = vz_client.elem(t_client_idxs);

 
    double DeltaVMinima = -1.0; 

    int sols_tot = 0; 

    // for every possible departure time from the service depot
    for (int t_ser=0 ; t_ser < available_t_service.n_elem ; t_ser++){



        // fetch valid arrival times and calculate possible times of flight
        arma::uvec valid_arrivals = arma::find(available_t_client > available_t_service(t_ser));
        arma::vec possible_tofs = available_t_client.elem(valid_arrivals) - available_t_service(t_ser);
        
        // service station velocities at departure time
        arma::vec3 v_service_loop = {available_vx_service(t_ser), available_vy_service(t_ser), available_vz_service(t_ser)};

        // r vec at departure time 
        arma::vec3 r_service = {available_x_service(t_ser), available_y_service(t_ser), available_z_service(t_ser)}; 

        // possible arrival locations and velocities of the client
        arma::vec valid_arrivals_x_client = available_x_client.elem(valid_arrivals);
        arma::vec valid_arrivals_y_client = available_y_client.elem(valid_arrivals);
        arma::vec valid_arrivals_z_client = available_z_client.elem(valid_arrivals);

        arma::vec valid_vx_client = available_vx_client.elem(valid_arrivals);
        arma::vec valid_vy_client = available_vy_client.elem(valid_arrivals);
        arma::vec valid_vz_client = available_vz_client.elem(valid_arrivals);

        // for all possible times of flight for the given departure time
        for(int t_cl = 0; t_cl < possible_tofs.n_elem ; t_cl++){
            // we don't need to check all possible times 

            // get the position vector of the client at time of arrival
            arma::vec3 r_client_loop = {valid_arrivals_x_client(t_cl),valid_arrivals_y_client(t_cl),valid_arrivals_z_client(t_cl)}; 
            arma::vec3 v_client_loop = {valid_vx_client(t_cl),valid_vy_client(t_cl),valid_vz_client(t_cl)};

            // calculate the lambert transfers that are possible
            std::vector<arma::vec> sols = lambert_solver(r_service, r_client_loop, possible_tofs(t_cl), MU_EARTH, 0, 1);
            
            // for all possible lambert transfers compute the delta V
            for (int sol=0 ; sol < sols.size(); sol++){
                sols_tot+=1;
                //std::cout << "Solution number: " << sols_tot << std::endl;

                arma::vec v1sol_loop = get_v1(sols,sol);
                arma::vec v2sol_loop = get_v2(sols,sol);

                //std::cout << " V1: " << get_v1(sols,sol);
                //std::cout << " V2: " << get_v1(sols,sol);

                arma::vec deltaV1 = v1sol_loop - v_service_loop;
                arma::vec deltaV2 = v2sol_loop - v_client_loop;

                double deltaTotal = arma::norm(deltaV1) + arma::norm(deltaV2);

                // to parallelise with omp you will have to remove this conditional.  
                // instead write solutions to a file or shared array and have a sparate loop identify the optimal solution 
                if ((sols_tot > 1) && (abs(deltaTotal) < DeltaVMinima) ){

                    DeltaVMinima = abs(deltaTotal);

                    v1sol = v1sol_loop;
                    v2sol = v2sol_loop; 

                    r1sol = r_service;
                    r2sol = r_client_loop;

                    tof_optimal = possible_tofs(t_cl);

                }

                else if (sols_tot == 1)
                {
                    DeltaVMinima = abs(deltaTotal);
                }
                


            }
            

        }


    }


}