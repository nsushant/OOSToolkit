
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include <chrono>
#include <ratio>
#include <thread>

#include "Nbody.hpp"
#include "LambertSolver.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include "Local_search.hpp"

int main(int argc, char *argv[]){
    double altitude_km = 700.0;
    int num_planes = 5;   
    int num_satellites = 50;  
    int relative_phase = 1;   

    double altitude = R_EARTH + altitude_km*1000.0;
    double period = 2*M_PI * sqrt(pow(altitude,3) / MU_EARTH);
    double dt = 10.0;
    double t_final = 60000;

    std::cout << "orbital period: " << period;  

    force_model force_options;
    force_options.includeJ2 = false;
    force_options.includeMutual = false;

    std::vector<satellite_object> sats = build_walker_constellation(num_planes, num_satellites, relative_phase,
                                                        altitude_km * 1000.0,

                                                        56.0 * M_PI/180.0,

                                                        100.0);

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
    std::ofstream csv("../data/WalkerDelta.csv");
    csv << "time_s,index,name,x,y,z,vx,vy,vz\n";


    double t = 0.0;
    int step = 0;


    std::cout<<"Running Simulation"; 
    while(t < (t_final - 1e-9)){
        // only write outputs to csv files every 10 timesteps
        // I am just thinning the datafile here so we can store it easily 

        if(step % 10 == 0){

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
        //showProgressBar(t, t_final);
        

    }
    
    std::cout << std::endl << "Done!" << std::endl;

    csv.close();
    std::cout << "CSV saved to data/WalkerDelta.csv\n";

    /*
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
    std::cout << std::setprecision(10) << "V1 (at service station) (m/s): " << get_v1(lambert_solutions,sols) << std::endl;
    std::cout << std::setprecision(10) << "V2 (at client) (m/s): " << get_v2(lambert_solutions,sols) << std::endl;

    }*/


    arma::vec v1sol,v2sol,r1sol,r2sol;
    std::vector<arma::vec> trajs; 

    double tof_optimal; 
    
    DataFrame simfile("../data/WalkerDelta.csv");    

    //std::vector<double> t_final_trys = {100,500,1000,1500,10000,15000,20000};  

    //for (int rep = 0; rep < t_final_trys.size(); rep++){

    //std::cout << "t_final = " << t_final_trys[rep]<<"-----------"<<std::endl; 
    
    //const auto start = std::chrono::high_resolution_clock::now(); 

    double DeltaVMinima; 

    find_optimal_trajectory("sat_10","service_1", 0.0, 60000.0, v1sol, v2sol, r1sol, r2sol,tof_optimal, trajs, simfile,"lambert",true,DeltaVMinima);    
    
    std::ofstream traj_file("../data/Trajectories.csv");
    //deltaTotal,possible_tofs(t_cl), available_t_service(t_ser),available_t_client(t_cl)
    traj_file << "deltaV,v1,v2,tof,t_depot,t_client,x_ser,y_ser,z_ser,x_cl,y_cl,z_cl,sol_num,vx_depot,vy_depot,vz_depot,vx_client,vy_client,vz_client\n";


    for(int it = 0 ; it < trajs.size() ; it++){

        arma::vec row = trajs[it]; 

        traj_file << row(0) << "," << row(1) << "," << row(2) << ","
                    << row(3) << "," << row(4) << "," << row(5) << ","
                    << row(6) << "," << row(7) << ","  << row(8) << "," 
                    << row(9) << "," << row(10) << ","<< row(11) << "," << row(12)
                    << "," << row(13) << "," << row(14)<< "," << row(15)
                    <<"," << row(16) << "," << row(17)<< "," << row(18)<<"\n";
        
    }

    traj_file.close();

    std::vector<std::string> satnames = simfile["name"];

    std::string depot_name = "service_1"; 
    std::string cl1_name = "sat_0"; 
    std::string cl2_name = "sat_3";
    std::string cl3_name = "sat_4";

    arma::vec t_sat = simfile.getNumeric("time_s"); 

    arma::uvec service_sat_idxs = find_idxs_of_match(satnames,depot_name);

    arma::uvec cl1 = find_idxs_of_match(satnames,cl1_name);
    arma::uvec cl2 = find_idxs_of_match(satnames,cl2_name);
    arma::uvec cl3 = find_idxs_of_match(satnames,cl3_name);


    
    arma::vec t_depot = t_sat.elem(service_sat_idxs);  
    arma::vec t_cl1 = t_sat.elem(cl1);  
    arma::vec t_cl2 = t_sat.elem(cl2);  
    arma::vec t_cl3 = t_sat.elem(cl3);  

    // initialization of schedule 
    std::vector<double> t_depart = {0, 8000, 16000, 24000,34000,42000,0.0}; 
    std::vector<double> t_arrive = {0, 6000,14000,22000,32000,40000,50000};

    double deltaV_of_schedule_init;
    
    
    /*
    
    // This part would go through all possible arrival and departure times
  
    std::vector<std::string> sat_names_in_schedule = {depot_name,cl1_name,depot_name,cl2_name,depot_name,"sat_10",depot_name};     
    
    
    const auto start = std::chrono::high_resolution_clock::now();
    schedule_struct init_schedule = create_schedule(deltaV_of_schedule_init, t_arrive, t_depart, sat_names_in_schedule,simfile);
    
    std::cout<<"Initial Schedule"<<std::endl;
    view_schedule(init_schedule); 
    
    schedule_struct findopt_schedule = local_search_opt_schedule(deltaV_of_schedule_init, init_schedule, 100, simfile);
    const auto stop = std::chrono::high_resolution_clock::now();
    auto calculation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout<<"Optimal Schedule"<<std::endl; 
    view_schedule(findopt_schedule); 
    std::cout << "execution time (ms): " << calculation_duration.count() << "\n";

    */
    
    

    std::vector<std::string> sat_names_in_schedule = {depot_name,cl1_name,depot_name,cl2_name,depot_name,"sat_10",depot_name};     

    
    const auto start = std::chrono::high_resolution_clock::now();

    double service_time = 1000;
    schedule_struct init_schedule = create_schedule_lambert_only(deltaV_of_schedule_init,t_arrive, t_depart, sat_names_in_schedule,simfile, service_time);

    //double total_t = (init_schedule.blocks[-1].arrival_time - init_schedule.blocks[0].arrival_time);  

    std::cout<<" Initial Schedule "<<std::endl;
    std::cout<<" Total Delta V: "<< deltaV_of_schedule_init << std::endl;
    view_schedule(init_schedule); 
    std::cout<<"\n";


    //std::vector<std::string> moves_for_local_search = {"add arrival", "sub arrival", "add departure", "sub departure"};
    
    std::vector<std::string> moves_for_local_search = {"move_dt2"};

    
    std::vector<double> grid_move_dts = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800}; 

    std::vector<double> deltaVsobtained; 

    std::vector<double> TimeImprovementObtained; 


    for (int dts = 0 ; dts < grid_move_dts.size(); dts++) {

    double deltaV_of_schedule_init2 = deltaV_of_schedule_init;
    
    const auto start = std::chrono::high_resolution_clock::now();
    schedule_struct findopt_schedule = local_search_opt_schedule_lambert_only(deltaV_of_schedule_init2, init_schedule, grid_move_dts[dts],simfile,service_time,moves_for_local_search);
    const auto stop = std::chrono::high_resolution_clock::now();
    
    auto calculation_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "execution time (ms): " << calculation_duration.count() << "\n";

    double arrival_improved = init_schedule.blocks[init_schedule.blocks.size()-1].arrival_time- findopt_schedule.blocks[findopt_schedule.blocks.size()-1].arrival_time; 
    
    double deltaV_improved = deltaV_of_schedule_init - deltaV_of_schedule_init2; 

    view_schedule(findopt_schedule);

    std::cout<<"Total Delta V: "<< deltaV_of_schedule_init2 << std::endl;

    TimeImprovementObtained.push_back(arrival_improved);
    deltaVsobtained.push_back(deltaV_improved);

    }
    

    std::ofstream DeltaVsGrid("../data/DeltaV_vs_movesize.csv");
    DeltaVsGrid << "move_size,deltaV_improvement,time_improvement"<<"\n";  

    //write rows to csv
 
    for(int row= 0 ; row < grid_move_dts.size() ; row++ ){

        DeltaVsGrid<<grid_move_dts[row]<<","<<deltaVsobtained[row]<<","<<TimeImprovementObtained[row]<<"\n";
    
    }
    

    DeltaVsGrid.close(); 

    


    return 0;


    }       


    












