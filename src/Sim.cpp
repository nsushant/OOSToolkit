
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

  run_simulation("WalkerDelta.csv", "walker_delta", 550000, 10,
                 altitude_m, num_planes, num_satellites, relative_phase,
                 deg_to_rads(inclination_in_deg), 1.0, deg_to_rads(inclinationDiff_in_deg));

  std::cout << "finished running sim" << "\n";

  std::cout << "------------------Running Local Search---------------------" << std::endl; 

  // Setting up an initial schedule
  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  // inputs to generate initial schedule
  int num_sat_visits = 3;
  double t_final = 50000;
  double service_time = 1000;

  std::string depot_name = "service_1";
  std::vector<std::string> client_satnames = {"sat_0", "sat_3", "sat_10", "sat_1", "sat_4", "sat_12", "sat_15", "sat_20", "sat_14", "sat_18"};

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
  std::vector<std::string> moves_to_consider = {"add departure", "sub arrival", "sub departure", "add arrival", "move_sub_traj", "move_add_traj"}; //, "move_dt2", "move_dt2_inv"};

  // run_local_search_tfixed

  schedule_struct ls = run_local_search(  simfile, move_size, moves_to_consider,
                                          sat_names_in_schedule, t_depart, t_arrive,
                                          deltaV_of_schedule_heuristic, service_time);

  view_schedule(ls);

  // vn_search(deltaV_of_schedule_init, ls, move_size,simfile, service_time, moves_to_consider, 500);

  // run_vn_search_fixed_tarrive(simfile, move_size, moves_to_consider,
  //             sat_names_in_schedule, t_depart, t_arrive,
  //           deltaV_of_schedule_heuristic, service_time, 400);

  std::cout << "Heuristic Delta V : " << deltaV_of_schedule_heuristic << "\n";

  std::cout << "------------------Running Exact Method---------------------" << std::endl;

  /*
  double initdeltavDP;
  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  std::vector<int> num_sat_visits_arr = {1, 2, 3, 4, 5, 6, 7, 8};

  std::vector<std::chrono::duration<double, std::milli>> scalingts;

  for (int v : num_sat_visits_arr)
  {
    int num_sat_visits = v;
    double t_final = 150000 * v / 8 - (150000 * v / 8 % 100);
    double service_time = 1000;

    std::string depot_name = "service_1";
    std::vector<std::string> client_satnames = {"sat_0", "sat_3", "sat_10", "sat_1", "sat_4", "sat_12", "sat_15", "sat_20", "sat_14", "sat_18"};

    // initializing init schedule
    std::vector<double> t_depart;
    std::vector<double> t_arrive;

    init_dep_arrival_times_strict_timespan(t_depart, t_arrive, t_final, service_time, num_sat_visits);
    std::vector<std::string> sat_names_in_schedule = init_satname_array("service_1", client_satnames, false, num_sat_visits);

    schedule_struct schedule_init = create_schedule_lambert_only(initdeltavDP, t_arrive, t_depart, sat_names_in_schedule, simfile, service_time);

    auto start = std::chrono::high_resolution_clock::now();

    double deltaVoptimal_ex_search = run_es(schedule_init);

    auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "\n";
    std::cout << "Num Visits: " << v << "\n";
    std::cout << "------milliseconds : " << duration.count() << "\n";

    scalingts.push_back(duration);

  }

  */
  double initdeltavDP;

  schedule_struct schedule_init = create_schedule_lambert_only(initdeltavDP, t_arrive, t_depart, sat_names_in_schedule, simfile, service_time);

  double deltaVoptimal_ex_search = run_es(schedule_init);

  std::cout << "Exact method DeltaV: " << deltaVoptimal_ex_search << "\n";
  std::cout << "Gap : " << std::abs(deltaV_of_schedule_heuristic - deltaVoptimal_ex_search) / deltaVoptimal_ex_search * 100 << " % " << "\n";

  // int ret = dynamic_program_fixed_tasksize_Tfixed(schedule_init, 100, 1,schedule_init.blocks.size()-1, 70000,simfile);

  // view_schedule(schedule_init);

  /*
  std::cout << (int)schedule_init.blocks.size() << std::endl;

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
