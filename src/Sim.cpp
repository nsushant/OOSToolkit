
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
#include "VNS.cpp"





void run_es(schedule_struct sched_in){


    for(int i=1 ; i < sched_in.blocks.size(); i++){
        
        task_block fromblock = sched_in.blocks[i-1]; 
        task_block toblock = sched_in.blocks[i];

        run_exhaustive_search(fromblock.satname, toblock.satname, fromblock.departure_time, toblock.arrival_time, "WalkerDelta.csv", "","lambert");
        std::cout<<"\n"; 

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

  run_simulation( "WalkerDelta.csv", "walker_delta", 100000,  10,
                 altitude_m, num_planes, num_satellites, relative_phase,
               deg_to_rads(inclination_in_deg), 1.0,deg_to_rads(inclinationDiff_in_deg));



  std::cout<<"finished running sim"<<"\n";


  //std::vector<double> move_size;

  //for (int i = 0; i<=6000 ; i+=100){
    
     // move_size.push_back(i);
    
  //}

  

  std::cout << "------------------Running Local Search---------------------" << std::endl;

  DataFrame simfile("../data/WalkerDelta.csv");

  std::vector<std::string> satnames = simfile["name"];

  // inputs to generate initial schedule
  int num_sat_visits = 3;
  double t_final = 100000;
  double service_time = 1000;
  std::string depot_name = "service_1";
  std::vector<std::string> client_satnames = {"sat_0", "sat_3", "sat_10","sat_1","sat_4","sat_12","sat_15","sat_20","sat_14","sat_18"};

  // initializing init schedule
  std::vector<double> t_depart;
  std::vector<double> t_arrive;

  init_dep_arrival_times_strict_timespan(t_depart, t_arrive, t_final, service_time, num_sat_visits);
  std::vector<std::string> sat_names_in_schedule = init_satname_array("service_1", client_satnames, false, num_sat_visits);

  double deltaV_of_schedule_init;

  // ideal case
  std::vector<double> move_size = {100, 200, 400,500,600, 700, 1000, 1500};

  //std::vector<double> move_size = {100,1000};

  arma::vec x_sat = simfile.getNumeric("x");
  arma::vec y_sat = simfile.getNumeric("y");
  arma::vec z_sat = simfile.getNumeric("z");
  arma::vec vx_sat = simfile.getNumeric("vx");
  arma::vec vy_sat = simfile.getNumeric("vy");
  arma::vec vz_sat = simfile.getNumeric("vz");

  std::vector<std::string> sats_in_demand = {"sat_0", "sat_3", "sat_10", "service_1"};

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

  double cooling_param = 0.99;
  int maxiter = 80;
  double T_init = 90;
  // formultaion 1
  // std::vector<std::string> moves_to_consider = {"sub arrival"};

  // formulation 2
  std::vector<std::string> moves_to_consider = {"add departure", "sub arrival", "sub departure", "add arrival","move_dt2"};

  // formulation 3
  // std::vector<std::string> moves_for_local_search = {"swap slots"};

  schedule_struct ls = run_local_search( simfile, move_size,moves_to_consider,
                                        sat_names_in_schedule, t_depart, t_arrive,
                                       deltaV_of_schedule_init, service_time );

    
  view_schedule(ls); 
  //std::cout << "running tabu search" << std::endl; 

  //run_tabu_search(simfile, move_size,moves_to_consider,
    //                sat_names_in_schedule,t_depart, t_arrive,
      //              deltaV_of_schedule_init, service_time);
 
  
  std::cout << "running variable neighbourhood search" << std::endl; 

  //vn_search(deltaV_of_schedule_init, ls, move_size,simfile, service_time, moves_to_consider, 500); 

  //view_schedule(ls); 

  // run_simulated_annealing( simfile, move_size,
  //                         moves_to_consider,
  //                       sat_names_in_schedule ,
  //                     t_depart, t_arrive,
  //                  deltaV_of_schedule_init, service_time, cooling_param, maxiter, T_init);

  //run_vn_search(simfile, move_size, moves_to_consider,
    //            sat_names_in_schedule, t_depart, t_arrive,
      //          deltaV_of_schedule_init, service_time,1000);
 
  std::cout << "------------------Running Dynamic Program---------------------" << std::endl;

  double initdeltavDP;


  // schedule_struct schedule_init =create_schedule(initdeltavDP,t_arrive,t_depart,sat_names_in_schedule,simfile) ;

  schedule_struct schedule_init = create_schedule_lambert_only(initdeltavDP, t_arrive, t_depart, sat_names_in_schedule, simfile, service_time);
    

  run_es(schedule_init); 

  view_schedule(schedule_init); 

  //int ret = dynamic_program_fixed_tasksize_Tfixed(schedule_init, 100, 1,schedule_init.blocks.size()-1, 70000,simfile);


  //view_schedule(schedule_init);

  /*  
  std::cout << (int)schedule_init.blocks.size() << std::endl;

  finding_individual_minimas_dynamic_programming(schedule_init, simfile, 100);

  view_schedule(schedule_init);

  double deltaVoptimal_exact = 0.0;

  for (int b = 0; b < schedule_init.blocks.size(); b++)
  {

    deltaVoptimal_exact += schedule_init.blocks[b].deltaV_arrival;
  }

  std::cout << "Exact method DeltaV: " << deltaVoptimal_exact << "\n";
  std::cout << "Gap : " << std::abs(deltaV_of_schedule_init - deltaVoptimal_exact) / deltaVoptimal_exact * 100 << " %" << "\n";


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
