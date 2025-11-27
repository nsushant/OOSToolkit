
#include <armadillo>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ratio>
#include <string>
#include <thread>
#include <vector>

#include "LambertSolver.hpp"
#include "Local_search.hpp"
#include "Nbody.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include "Exact_methods.hpp"
#include "Simulated_Annealing.hpp"

int main(int argc, char *argv[]) {


  std::cout << "------------------Running Simulation---------------------"<<std::endl; 

  double altitude_m = 700*1000; 
  int num_planes = 5;
  int num_satellites = 30;
  int relative_phase = 1;
  double inclination_in_deg = 56; 
  double inclinationDiff_in_deg = 2;    

  //run_simulation( "WalkerDelta.csv", "walker_delta", 60000,  10, 
    //              altitude_m, num_planes, num_satellites, relative_phase,
      //            deg_to_rads(inclination_in_deg), 1.0,deg_to_rads(inclinationDiff_in_deg));
  

  //std::cout << "------------------Running exact method (DP)---------------------"<<std::endl; 


  

  //run_exhaustive_search("service_1", "sat_0", 0.0, 30000.0, "WalkerDelta.csv", "Trajectories.csv","lambert"); 

  //run_exhaustive_search("service_1", "sat_0", 0.0, 30000.0, "WalkerDelta.csv", "","edelbaum"); 

  
  std::cout << "------------------Running Local Search---------------------"<<std::endl; 
  

  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];
  
  // inputs to generate initial schedule
  int num_sat_visits = 3; 
  double t_final = 50000; 
  double service_time = 1000;
  std::string depot_name = "service_1"; 
  std::vector<std::string> client_satnames  = {"sat_0", "sat_3", "sat_10"};
  
  // initializing init schedule 
  std::vector<double> t_depart; 
  std::vector<double> t_arrive; 

  init_dep_arrival_times_strict_timespan(t_depart, t_arrive, t_final, service_time, num_sat_visits);
  std::vector<std::string> sat_names_in_schedule = init_satname_array("service_1", client_satnames, false, num_sat_visits);

  double deltaV_of_schedule_init;

  std::vector<double> move_size = {100,200,500,1000,1500,3000}; 

  arma::vec x_sat = simfile.getNumeric("x");
  arma::vec y_sat = simfile.getNumeric("y");
  arma::vec z_sat = simfile.getNumeric("z");
  arma::vec vx_sat = simfile.getNumeric("vx");
  arma::vec vy_sat = simfile.getNumeric("vy");
  arma::vec vz_sat = simfile.getNumeric("vz");


  std::vector<std::string> sats_in_demand= {"sat_0" , "sat_3" , "sat_10" , "service_1"}; 

  for(std::string sname : sats_in_demand){

    arma::uvec idxs_sats = find_idxs_of_match(satnames,sname);
    
    arma::vec x = x_sat.elem(idxs_sats);
    arma::vec y = y_sat.elem(idxs_sats);
    arma::vec z = z_sat.elem(idxs_sats);


    arma::vec vx = vx_sat.elem(idxs_sats);
    arma::vec vy = vy_sat.elem(idxs_sats);
    arma::vec vz = vz_sat.elem(idxs_sats);
    
    arma::vec r  = { x(0), y(0),z(0)};
    arma::vec v  = {vx(0) , vy(0) , vz(0)};

    
    orbital_elements eleminit = orb_elems_from_rv(r,v);


    std::cout<<"sat: "<<sname<<std::endl;
    std::cout<<"inclination: "<<eleminit.inclination<<std::endl; 
    std::cout<<"RAAN: "<<eleminit.RAAN<<std::endl;
    std::cout<<"SemiMajor Axis: "<< eleminit.semi_major_axis<<std::endl;
    std::cout<<"__________________________________________"<<std::endl;

  }
    

  double cooling_param = 0.99; 
  int maxiter = 40;
  double T_init = 50; 
  //formultaion 1 
  //std::vector<std::string> moves_to_consider = {"sub arrival"};


  // formulation 2 
  std::vector<std::string> moves_to_consider = { "add departure","sub arrival", "add arrival", "sub departure"};


  // formulation 3  
  //std::vector<std::string> moves_for_local_search = {"swap slots"};

  run_local_search( simfile, move_size,moves_to_consider,
                    sat_names_in_schedule, t_depart, t_arrive,
                    deltaV_of_schedule_init, service_time );

  /*
  run_simulated_annealing( simfile, move_size,  
                            moves_to_consider,
                            sat_names_in_schedule ,
                            t_depart, t_arrive, 
                           deltaV_of_schedule_init, service_time, cooling_param, maxiter, T_init);
  */

  std::cout << "------------------Running Dynamic Program---------------------"<<std::endl; 

  double initdeltavDP; 
  //schedule_struct schedule_init =create_schedule(initdeltavDP,t_arrive,t_depart,sat_names_in_schedule,simfile) ;
  
  schedule_struct schedule_init =create_schedule_lambert_only(initdeltavDP,t_arrive,t_depart, sat_names_in_schedule,simfile, service_time);

  std::cout<<(int)schedule_init.blocks.size() << std::endl;
  
  finding_individual_minimas_dynamic_programming(schedule_init, simfile, 100);

  view_schedule(schedule_init);


  double deltaVoptimal_exact = 0.0; 

  for(int b = 0 ; b < schedule_init.blocks.size() ; b++){

    deltaVoptimal_exact += schedule_init.blocks[b].deltaV_arrival; 


  }   

  std::cout<<"Exact method DeltaV: " << deltaVoptimal_exact << "\n";
  std::cout<<"Gap : "<< std::abs(deltaV_of_schedule_init - deltaVoptimal_exact)/deltaVoptimal_exact * 100 <<" %" << "\n"; 

  /*

  std::ofstream DeltaVsGrid("../data/DeltaV_vs_movesize.csv");
  DeltaVsGrid << "move_size,deltaV_improvement,time_improvement" << "\n";

  // write rows to csv

  for (int row = 0; row < grid_move_dts.size(); row++) {

    DeltaVsGrid << grid_move_dts[row] << "," << deltaVsobtained[row] << ","
                << TimeImprovementObtained[row] << "\n";
  }

  DeltaVsGrid.close();
  */

  return 0;


}
