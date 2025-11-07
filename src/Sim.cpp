
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
  double inclination_in_deg = 56; 

  run_simulation( "WalkerDelta.csv", "walker_delta", 60000,  10, 
                  altitude_m, num_planes, num_satellites, relative_phase,
                  deg_to_rads(inclination_in_deg), 1.0);
  

  std::cout << "------------------Running exhaustive search---------------------"<<std::endl; 


  run_exhaustive_search("service_1", "sat_0", 0.0, 30000.0, "WalkerDelta.csv", "Trajectories.csv"); 

  
  std::cout << "------------------Running Local Search---------------------"<<std::endl; 
  
  
  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  std::string depot_name = "service_1";
  std::string cl1_name = "sat_0";
  std::string cl2_name = "sat_3";
  std::string cl3_name = "sat_4";

  arma::vec t_sat = simfile.getNumeric("time_s");

  arma::uvec service_sat_idxs = find_idxs_of_match(satnames, depot_name);

  arma::uvec cl1 = find_idxs_of_match(satnames, cl1_name);
  arma::uvec cl2 = find_idxs_of_match(satnames, cl2_name);
  arma::uvec cl3 = find_idxs_of_match(satnames, cl3_name);

  arma::vec t_depot = t_sat.elem(service_sat_idxs);
  arma::vec t_cl1 = t_sat.elem(cl1);
  arma::vec t_cl2 = t_sat.elem(cl2);
  arma::vec t_cl3 = t_sat.elem(cl3);

  // initialization of schedule
  std::vector<double> t_depart = {0, 8000, 16000, 24000,34000,42000,0.0}; 
  std::vector<double> t_arrive = {0, 6000,14000,22000,32000,40000,50000};
  //std::vector<double> deadline_array = {0, 32000, 0, 40000, 0, 50000, 0};

  

  double deltaV_of_schedule_init;

  /*


  */
  std::vector<std::string> sat_names_in_schedule = {depot_name,cl1_name,depot_name,cl2_name,depot_name,"sat_10",depot_name};     

  //std::vector<std::string> sat_names_in_schedule = {
   //   depot_name, "sat_0",  depot_name, "sat_1",  depot_name, "sat_3",
     // depot_name, "sat_5",  depot_name, "sat_7",  depot_name, "sat_9",
      //depot_name, "sat_10", depot_name, "sat_20", depot_name};


  double service_time = 1000;
  schedule_struct init_schedule = create_schedule_lambert_only(
      deltaV_of_schedule_init, t_arrive, t_depart, sat_names_in_schedule,
      simfile, service_time);


  std::cout << "init schedule"<<"\n";
  // double total_t = (init_schedule.blocks[-1].arrival_time -
  // init_schedule.blocks[0].arrival_time);

  std::cout << " Initial Schedule " << std::endl;
  std::cout << " Total Delta V: " << deltaV_of_schedule_init << std::endl;
  view_schedule(init_schedule);
  std::cout << "\n";

  // std::vector<std::string> moves_for_local_search = {"add arrival", "sub
  // arrival", "add departure", "sub departure"};
  
  // formulation 1   
  std::vector<std::string> moves_for_local_search = {"sub arrival", "add departure"};

  // formulation 2 
  //std::vector<std::string> moves_for_local_search = {"sub arrival", "add departure"};


  // formulation 3  


  std::vector<double> grid_move_dts = {
      600}; // 100,  200,  300,  400,  500,  600,
             // 700,  800,  900,  1000, 1100, 1200,
             // 1300, 1400, 1500, 1600, 1700, 1800};

  std::vector<double> deltaVsobtained;

  std::vector<double> TimeImprovementObtained;

  for (int dts = 0; dts < grid_move_dts.size(); dts++) {

    double deltaV_of_schedule_init2 = deltaV_of_schedule_init;

    const auto start = std::chrono::high_resolution_clock::now();
    schedule_struct findopt_schedule = local_search_opt_schedule_lambert_only(
        deltaV_of_schedule_init2, init_schedule, grid_move_dts[dts], simfile,
        service_time, moves_for_local_search);
    const auto stop = std::chrono::high_resolution_clock::now();

    auto calculation_duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "execution time (ms): " << calculation_duration.count()
              << "\n";

    double arrival_improved =
        init_schedule.blocks[init_schedule.blocks.size() - 1].arrival_time -
        findopt_schedule.blocks[findopt_schedule.blocks.size() - 1]
            .arrival_time;

    double deltaV_improved = deltaV_of_schedule_init - deltaV_of_schedule_init2;

    view_schedule(findopt_schedule);

    std::cout << "Total Delta V: " << deltaV_of_schedule_init2 << std::endl;

    TimeImprovementObtained.push_back(arrival_improved);
    deltaVsobtained.push_back(deltaV_improved);
  }

  std::ofstream DeltaVsGrid("../data/DeltaV_vs_movesize.csv");
  DeltaVsGrid << "move_size,deltaV_improvement,time_improvement" << "\n";

  // write rows to csv

  for (int row = 0; row < grid_move_dts.size(); row++) {

    DeltaVsGrid << grid_move_dts[row] << "," << deltaVsobtained[row] << ","
                << TimeImprovementObtained[row] << "\n";
  }

  DeltaVsGrid.close();

  return 0;
}
