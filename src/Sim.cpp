
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

int main(int argc, char *argv[]) {


  std::cout << "------------------Running Simulation---------------------"<<std::endl; 

  double altitude_m = 700*1000; 
  int num_planes = 5;
  int num_satellites = 20;
  int relative_phase = 1;
  double inclination_in_deg = 1; //56 
  double inclinationDiff_in_deg = 2;    

  run_simulation( "WalkerDelta.csv", "circular_orbits", 60000,  10, 
                  altitude_m, num_planes, num_satellites, relative_phase,
                  deg_to_rads(inclination_in_deg), 1.0,deg_to_rads(inclinationDiff_in_deg));
  

  std::cout << "------------------Running exhaustive search---------------------"<<std::endl; 


  //run_exhaustive_search("service_1", "sat_0", 0.0, 30000.0, "WalkerDelta.csv", "Trajectories.csv","edelbaum"); 

  run_exhaustive_search("service_1", "sat_0", 0.0, 30000.0, "WalkerDelta.csv", "","edelbaum"); 

  
  std::cout << "------------------Running Local Search---------------------"<<std::endl; 
  

  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  std::vector<double> t_depart = {0, 8000, 16000, 24000,34000,42000,0.0}; 
  std::vector<double> t_arrive = {0, 6000,14000,22000,32000,40000,50000};

  double deltaV_of_schedule_init;


  std::string depot_name = "service_1"; 
  

  std::vector<std::string> sat_names_in_schedule = {depot_name,"sat_0",depot_name,"sat_3",depot_name,"sat_10",depot_name};     

  double service_time = 1000;
  double move_size = 600; 

  //formultaion 1 
  std::vector<std::string> moves_to_consider = {"sub arrival"};


  // formulation 2 
  //std::vector<std::string> moves_for_local_search = {"sub arrival", "add departure"};


  // formulation 3  
  //std::vector<std::string> moves_for_local_search = {"swap slots"};



  run_local_searh(  simfile, move_size,moves_to_consider,
                    sat_names_in_schedule, t_depart, t_arrive, 
                    deltaV_of_schedule_init, service_time); 


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
