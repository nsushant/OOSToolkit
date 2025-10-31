
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
#include "Nbody.hpp"
#include "Trajectory_selection.hpp"


/*
Author : S. Nigudkar (2025)


The aim of the below function is to choose lambert parameters r1,r2,tof and compute the associated delta Vs 
for the most efficient orbital transfers given a requested service deadline t_service, the start time
t_prev, the name of the client satellite (departure destination) and that of the service satellite (arrival destination).

*/



void find_optimal_trajectory(std::string service_satname,std::string client_satname, double t_prev, double t_request ,arma::vec &v1sol, arma::vec &v2sol, arma::vec &r1sol, arma::vec &r2sol, double &tof_optimal, std::vector<arma::vec>&trajs
,DataFrame simfile, std::string method, bool write_to_file, double &DeltaVMinima){

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

    DeltaVMinima = -1.0; 

    int sols_tot = 0; 
    
    //std::cout<<"finding optimal trajectory"; 
    // for every possible departure time from the service depot
    for ( int t_ser=0 ; t_ser < available_t_service.n_elem ; t_ser++ ){

        // fetch valid arrival times and calculate possible times of flight
        arma::uvec valid_arrivals = arma::find(available_t_client > available_t_service(t_ser));
        arma::vec possible_tofs = available_t_client.elem(valid_arrivals) - available_t_service(t_ser);
        arma::vec valid_t_client = available_t_client.elem(valid_arrivals);
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

            if (method == "lambert") {

            std::vector<arma::vec> sols = lambert_solver(r_service, r_client_loop, possible_tofs(t_cl), MU_EARTH, 0, 4);
            
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


                    if (write_to_file == true){
                    
                    arma::vec costs = {deltaTotal,arma::norm(v1sol_loop),arma::norm(v2sol_loop),possible_tofs(t_cl), available_t_service(t_ser),valid_t_client(t_cl)};
                    arma::vec sol_v = {(double)sol};
                    arma::vec sol_v2 = {v_service_loop(0),v_service_loop(1),v_service_loop(2),v_client_loop(0),v_client_loop(1),v_client_loop(2)};
                    arma::vec vec_append = arma::join_cols(costs,r_service); 
                    vec_append = arma::join_cols(vec_append,r_client_loop);
                    vec_append = arma::join_cols(vec_append,sol_v);
                    vec_append = arma::join_cols(vec_append,sol_v2);
                    trajs.push_back(vec_append);

                    }

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


            else if (method == "edelbaum")
            {
                double v_depot_edelbaum = arma::norm(v_service_loop);
                double v_client_edelbaum = arma::norm(v_client_loop);
                double plane_diff_angle; 
                
                double deltaV = calculate_edelbaum_deltaV(v_depot_edelbaum,v_client_edelbaum, plane_diff_angle);

        
            }
            
            

        }

        //showProgressBar(t_ser, (int)available_t_service.n_elem);
        //std::cout<<"\n";
    }


    //std::cout << "DeltaVMinima : " << DeltaVMinima << std::endl; 

}


double compact_optimal_calc(std::string satname1,std::string satname2,double departure_time, double arrival_time,DataFrame simfile){


    double deltaV_arrival, tof_optimal; 

    arma::vec v1sol, v2sol, r1sol, r2sol; 

    std::vector<arma::vec> trajs; 


    find_optimal_trajectory( satname1, satname2, departure_time, arrival_time, v1sol, v2sol, r1sol, r2sol,tof_optimal, trajs, simfile,"lambert",false,deltaV_arrival);

    return deltaV_arrival;


}



void find_optimal_trajectory_no_iter(std::string service_satname, std::string client_satname, double t_departure, double t_arrival, DataFrame simfile, double &DeltaVMinima){

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

    // Finding physical state params of the service satellite 
    arma::vec t_service_sat = t_sat.elem(service_sat_idxs);
    //arma::uword t_service_departure_idx = arma::index_min(arma::abs(t_service_sat - t_departure));

    arma::uvec idxs_greater_dep = arma::find(t_service_sat >= t_departure); 

    t_service_sat = t_service_sat.elem(idxs_greater_dep); 


    arma::uword t_service_departure_idx = arma::index_min(t_service_sat);


    double t_service_departure = t_service_sat(t_service_departure_idx);

    arma::vec x_service = x_sat.elem(service_sat_idxs);
    x_service = x_service.elem(idxs_greater_dep);
    double x_service_sat = x_service(t_service_departure_idx); 

    arma::vec y_service = y_sat.elem(service_sat_idxs);
    y_service = y_service.elem(idxs_greater_dep);
    double y_service_sat = y_service(t_service_departure_idx); 

    arma::vec z_service = z_sat.elem(service_sat_idxs);
    z_service = z_service.elem(idxs_greater_dep);
    double z_service_sat = z_service(t_service_departure_idx); 

    arma::vec vx_service = vx_sat.elem(service_sat_idxs);
    vx_service = vx_service.elem(idxs_greater_dep);
    double vx_service_sat = vx_service(t_service_departure_idx); 

    arma::vec vy_service = vy_sat.elem(service_sat_idxs);
    vy_service = vy_service.elem(idxs_greater_dep);
    double vy_service_sat = vy_service(t_service_departure_idx); 

    arma::vec vz_service = vz_sat.elem(service_sat_idxs);
    vz_service = vz_service.elem(idxs_greater_dep);
    double vz_service_sat = vz_service(t_service_departure_idx); 


    // Finding physical state params of the client satellite 
    arma::vec t_client = t_sat.elem(client_sat_idxs) ; 

    //arma::uword t_client_arrival_idx = arma::index_min(arma::abs(t_client - t_arrival));
    arma::uvec idxs_greater_arr = arma::find(t_client >= t_arrival); 
    t_client = t_client.elem(idxs_greater_arr); 

    arma::uword t_client_arrival_idx =  arma::index_min(t_client); 
    

    double t_client_arrival = t_client(t_client_arrival_idx); 

    arma::vec x_client = x_sat.elem(client_sat_idxs) ;
    x_client = x_client.elem(idxs_greater_arr);
    double x_client_sat = x_client(t_client_arrival_idx);    

    arma::vec y_client = y_sat.elem(client_sat_idxs) ;
    y_client = y_client.elem(idxs_greater_arr);
    double y_client_sat = y_client(t_client_arrival_idx);    

    arma::vec z_client = z_sat.elem(client_sat_idxs) ;
    z_client = z_client.elem(idxs_greater_arr);
    double z_client_sat = z_client(t_client_arrival_idx);    

    arma::vec vx_client = vx_sat.elem(client_sat_idxs) ;
    vx_client = vx_client.elem(idxs_greater_arr);
    double vx_client_sat = vx_client(t_client_arrival_idx);    

    arma::vec vy_client = vy_sat.elem(client_sat_idxs) ;
    vy_client = vy_client.elem(idxs_greater_arr);
    double vy_client_sat = vy_client(t_client_arrival_idx);    

    arma::vec vz_client = vz_sat.elem(client_sat_idxs) ;
    vz_client = vz_client.elem(idxs_greater_arr);
    double vz_client_sat = vz_client(t_client_arrival_idx);    

    DeltaVMinima = -1.0; 


    arma::vec r_service = {x_service_sat,y_service_sat,z_service_sat};
    arma::vec r_client = {x_client_sat,y_client_sat,z_client_sat};

    arma::vec v_service = {vx_service_sat,vy_service_sat,vz_service_sat}; 
    arma::vec v_client = {vx_client_sat,vy_client_sat,vz_client_sat}; 



    double tof_largest = t_client_arrival - t_service_departure; 
    std::vector<arma::vec> sols = lambert_solver(r_service, r_client, tof_largest, MU_EARTH, 0, 4);

    double deltaVsols = -1*pow(10,9); 

    
    for (int sol=0 ; sol < sols.size(); sol++){

        arma::vec v1sol_loop = get_v1(sols,sol);
        arma::vec v2sol_loop = get_v2(sols,sol);


        arma::vec deltaV1 = v1sol_loop - v_service;
        arma::vec deltaV2 = v2sol_loop - v_client;

        double deltaTotal = arma::norm(deltaV1) + arma::norm(deltaV2);


        if (deltaVsols > 0){

            if (deltaVsols > deltaTotal){
                
                deltaVsols = deltaTotal; 



            }


        }

        else{

            deltaVsols = deltaTotal; 

        }


    }

    DeltaVMinima = deltaVsols; 

}
