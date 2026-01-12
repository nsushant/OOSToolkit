
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

#include "Tabu_search.hpp"
#include "LambertSolver.hpp"
#include "Local_search.hpp"
#include "Nbody.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include "Exact_methods.hpp"
#include "Simulated_Annealing.hpp"
#include "VNS.hpp"

double run_es(schedule_struct &sched_in);

void run_comprehensive_scaling_tests();

double run_es(schedule_struct &sched_in)
{

  double DeltaVchange = 0.0;
  for (int i = 1; i < sched_in.blocks.size(); i++)
  {

    task_block fromblock = sched_in.blocks[i - 1];
    task_block toblock = sched_in.blocks[i];

    run_exhaustive_search(fromblock.satname, toblock.satname, fromblock.departure_time, toblock.arrival_constraint, DeltaVchange, "WalkerDelta.csv", "", "lambert");
    std::cout << "\n";
  }

  return DeltaVchange;
}

void show_orb_elems(DataFrame &simfile, std::vector<std::string> &sats_in_demand)
{

  std::vector<std::string> satnames = simfile["name"];

  arma::vec x_sat = simfile.getNumeric("x");
  arma::vec y_sat = simfile.getNumeric("y");
  arma::vec z_sat = simfile.getNumeric("z");
  arma::vec vx_sat = simfile.getNumeric("vx");
  arma::vec vy_sat = simfile.getNumeric("vy");
  arma::vec vz_sat = simfile.getNumeric("vz");

  for (std::string sname : sats_in_demand)
  {
    arma::uvec idxs_sats = find_idxs_of_match(satnames, sname);

    arma::vec x = x_sat.elem(idxs_sats);
    arma::vec y = y_sat.elem(idxs_sats);
    arma::vec z = z_sat.elem(idxs_sats);

    arma::vec vx = vx_sat.elem(idxs_sats);
    arma::vec vy = vy_sat.elem(idxs_sats);
    arma::vec vz = vz_sat.elem(idxs_sats);

    arma::vec r = {x(0), y(0), z(0)};
    arma::vec v = {vx(0), vy(0), vz(0)};

    orbital_elements eleminit = orb_elems_from_rv(r, v);

    std::cout << "sat: " << sname << std::endl;
    std::cout << "inclination: " << eleminit.inclination << std::endl;
    std::cout << "RAAN: " << eleminit.RAAN << std::endl;
    std::cout << "SemiMajor Axis: " << eleminit.semi_major_axis << std::endl;
    std::cout << "__________________________________________" << std::endl;
  }
}

int main(int argc, char *argv[])
{

  std::cout << "------------------Running Simulation---------------------" << std::endl;

  double altitude_m = 700 * 1000;
  int num_planes = 5;
  int num_satellites = 30;
  int relative_phase = 1;
  double inclination_in_deg = 56;
  double inclinationDiff_in_deg = 2;

  run_simulation("WalkerDelta.csv", "walker_delta", 770000, 10,
                 altitude_m, num_planes, num_satellites, relative_phase,
                 deg_to_rads(inclination_in_deg), 1.0, deg_to_rads(inclinationDiff_in_deg));

  std::cout << "finished running sim" << "\n";

  std::cout << "------------------Running Local Search---------------------" << std::endl; 

  // Setting up an initial schedule
  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  // inputs to generate initial schedule (reduced for testing)
  int num_sat_visits = 3;
  double t_final = 150000;
  double service_time = 1000;

  std::string depot_name = "service_1";
  std::vector<std::string> client_satnames = {"sat_0", "sat_3", "sat_10"};

  // initializing init schedule
  std::vector<double> t_depart;
  std::vector<double> t_arrive;

  init_dep_arrival_times_strict_timespan(t_depart, t_arrive, t_final, service_time, num_sat_visits);
  std::vector<std::string> sat_names_in_schedule = init_satname_array("service_1", client_satnames, false, num_sat_visits);

  show_orb_elems(simfile, sat_names_in_schedule);

  double deltaV_of_schedule_heuristic;

  // Running Local Search

  // ideal case (1 percent convergence)
  std::vector<double> move_size = {15000,10000,9000,8000,6000,5000, 4000, 3500, 3000, 2000, 1500, 1000, 500, 100}; //, 200, 400, 500, 600, 700, 1000, 1500, 2000, 2500, 2700, 3000};
  // formulation 1
  std::vector<std::string> moves_to_consider = {"add departure", "sub arrival", "sub departure", "add arrival", "move_dt2", "move_dt2_inv", "move_sub_traj", "move_add_traj"};

  // run_local_search_tfixed
  
  schedule_struct ls = run_local_search_tfixed(simfile, move_size, moves_to_consider,
                                          sat_names_in_schedule, t_depart, t_arrive,
                                          deltaV_of_schedule_heuristic, service_time);

  view_schedule(ls);
  
    
  // vn_search(deltaV_of_schedule_init, ls, move_size,simfile, service_time, moves_to_consider, 500);

  // run_vn_search_fixed_tarrive(simfile, move_size, moves_to_consider,
  //             sat_names_in_schedule, t_depart, t_arrive,
  //           deltaV_of_schedule_heuristic, service_time, 400);

  std::cout << "Heuristic Delta V : " << deltaV_of_schedule_heuristic << "\n";

  std::cout << "------------------Running Comprehensive Scaling Comparison---------------------" << std::endl;

  // Run the comprehensive scaling comparison test
  run_comprehensive_scaling_tests();

  // int ret = dynamic_program_fixed_tasksize_Tfixed(schedule_init, 100, 1,schedule_init.blocks.size()-1, 70000,simfile);

  // view_schedule(schedule_init);

  
  /*std::cout << (int)schedule_init.blocks.size() << std::endl;

  finding_individual_minimas_dynamic_programming(schedule_init, simfile, 100);

  view_schedule(schedule_init);

  double deltaVoptimal_exact = 0.0;

  for (int b = 0; b < schedule_init.blocks.size(); b++)
  {

    deltaVoptimal_exact += schedule_init.blocks[b].deltaV_arrival;
  }



  */

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

void run_comprehensive_scaling_tests()
{
  std::cout << "Starting comprehensive scaling comparison..." << std::endl;
  
  // Test configuration
  std::vector<int> visit_counts = {3, 4, 5, 6, 7, 8, 9, 10, 12};
  std::string depot_name = "service_1";
  std::vector<std::string> client_satnames = {"sat_0", "sat_3", "sat_10", "sat_1", "sat_4", "sat_12", "sat_15", "sat_20", "sat_14", "sat_18"};
  double service_time = 1000;
  
  // Setup CSV output
  std::ofstream results_file("../data/scaling_comparison_results.csv");
  results_file << "visits,method,time_ms,deltaV,iterations,quality_gap,success" << std::endl;
  
  DataFrame simfile("../data/WalkerDelta.csv");
  
  for (int visits : visit_counts)
  {
    std::cout << "\n===== Testing " << visits << " visits =====" << std::endl;
    
    // Setup instance
    double t_final = 150000 * visits / 8 - (150000 * visits / 8 % 100);
    std::vector<double> t_depart;
    std::vector<double> t_arrive;
    
    init_dep_arrival_times_strict_timespan(t_depart, t_arrive, t_final, service_time, visits);
    std::vector<std::string> sat_names_in_schedule = init_satname_array(depot_name, client_satnames, false, visits);
    
    double initdeltavDP;
    schedule_struct schedule_base = create_schedule_lambert_only(initdeltavDP, t_arrive, t_depart, sat_names_in_schedule, simfile, service_time);
    
    // Test Exact Method
    std::cout << "Testing Exact Method..." << std::endl;
    try {
      auto start_exact = std::chrono::high_resolution_clock::now();
      double deltaV_exact = run_es(schedule_base);
      auto end_exact = std::chrono::high_resolution_clock::now();
      
      double time_exact_ms = std::chrono::duration<double, std::milli>(end_exact - start_exact).count();
      
      std::cout << "Exact Method - Time: " << time_exact_ms << "ms, DeltaV: " << deltaV_exact << std::endl;
      
      // Write exact method results
      results_file << visits << ",exact," << time_exact_ms << "," << deltaV_exact << ",1,0.0,true" << std::endl;
      
      // Test Local Search (using same base schedule)
      std::cout << "Testing Local Search..." << std::endl;
      
      std::vector<double> move_size = {15000,10000,9000,8000,6000,5000, 4000, 3500, 3000, 2000, 1500, 1000, 500, 100};
      std::vector<std::string> moves_to_consider = {"add departure", "sub arrival", "sub departure", "add arrival", "move_dt2", "move_dt2_inv", "move_sub_traj", "move_add_traj"};
      
      auto start_ls = std::chrono::high_resolution_clock::now();
      double deltaV_ls = initdeltavDP;
      schedule_struct schedule_ls = run_local_search(simfile, move_size, moves_to_consider,
                                                   sat_names_in_schedule, t_depart, t_arrive,
                                                   deltaV_ls, service_time);
      auto end_ls = std::chrono::high_resolution_clock::now();
      
      double time_ls_ms = std::chrono::duration<double, std::milli>(end_ls - start_ls).count();
      
      // Calculate quality gap
      double quality_gap = std::abs(deltaV_ls - deltaV_exact) / deltaV_exact * 100.0;
      
      std::cout << "Local Search - Time: " << time_ls_ms << "ms, DeltaV: " << deltaV_ls 
                << ", Gap: " << quality_gap << "%" << std::endl;
      
      // Write local search results
      results_file << visits << ",local_search," << time_ls_ms << "," << deltaV_ls << ",0," 
                  << quality_gap << ",true" << std::endl;
                  
    } catch (const std::exception& e) {
      std::cout << "Error testing " << visits << " visits: " << e.what() << std::endl;
      results_file << visits << ",exact,0,0,0,0,false" << std::endl;
      results_file << visits << ",local_search,0,0,0,0,false" << std::endl;
    }
    
    std::cout << "Completed " << visits << " visits" << std::endl;
  }
  
  results_file.close();
  std::cout << "\nScaling comparison completed! Results saved to ../data/scaling_comparison_results.csv" << std::endl;
}
