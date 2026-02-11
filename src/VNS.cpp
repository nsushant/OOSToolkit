#include "Local_search.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include <armadillo>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <numeric>
#include <map>
#include <algorithm>
#include <limits>
/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/



int safe_sizet_to_int(size_t s)
{
  if (s > std::numeric_limits<int>::max())
    throw std::overflow_error("size_t too large for int");
  return static_cast<int>(s);
}

double find_minima_val(const std::vector<double>& v)
{
  if (v.empty()) {
    return std::numeric_limits<double>::infinity();
  }
  auto itmin = std::min_element(v.begin(), v.end());
  return *itmin; // smallest element from dereferenced iterator
}

int find_minima_index(const std::vector<double>& v)
{
  if (v.empty()) {
    return -1;
  }
  auto itmin = std::min_element(v.begin(), v.end());
  return std::distance(v.begin(), itmin);
}

void vn_search(double &init_deltaV, schedule_struct &init_schedule, std::vector<double> dt_move,
               DataFrame simfile, double service_time, std::vector<std::string>&move_methods, int max_iter)
{

  // Optimized feasibility check using departure_time <= arrival_constraint
  auto quick_feasibility_check = [&](const schedule_struct& schedule, int b) -> bool {
      // Service deadline constraint: departure_time <= arrival_constraint
      if (schedule.blocks[b].departure_time > schedule.blocks[b].arrival_constraint) {
        return false;
      }
      // Previous block continuity
      if (b > 0 && schedule.blocks[b].arrival_time <= schedule.blocks[b - 1].departure_time) {
        return false;
      }
      // Next block continuity
      if (b + 1 < schedule.blocks.size() && 
          schedule.blocks[b].departure_time >= schedule.blocks[b + 1].arrival_time) {
        return false;
      }
      // Positive time validation
      if (schedule.blocks[b].arrival_time <= 0 || schedule.blocks[b].departure_time <= 0) {
        return false;
      }
      return true;
  };

  int k = 0;
  int iterval = 0;

  double deltaVminima_so_far = init_deltaV;
  schedule_struct optimal_schedule;

  int move_choice = 0;

  std::random_device rd;          // used only once
  std::mt19937 gen(rd());         // good quality PRNG


  while (iterval < max_iter)
  {
    
    std::shuffle(move_methods.begin(), move_methods.end(), gen);

    iterval += 1;
    std::string move_current = move_methods[move_choice];
    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;
    std::vector<double> deltaVs_of_neighbourhood;

    double neighbourhood_minima;

    // loop over all moves
    for (int b = 1; b < init_schedule.blocks.size(); b++)
    {
      // Block index validation and move applicability check
      if (b >= init_schedule.blocks.size() || b < 0) {
        continue;
      }

      bool move_applicable = false;

      if ((move_current== "add" || move_current== "sub") && b > 0) {
        move_applicable = true;
      } else if ((move_current == "swap" || move_current == "swap_inv") && 
                 (b > 0 && b < init_schedule.blocks.size() - 1)) {
        move_applicable = true;
      }

      if (!move_applicable) {
        continue;
      }

      schedule_struct schedule_sol = init_schedule;

      // Validate dt values
      std::vector<double> dt = {0.5,0.3,0.2}; 
      
      std::srand(std::time(nullptr));
      
      int randomIndex = std::rand() % dt.size();

      double largest_step_size = init_schedule.blocks[b].arrival_time - 
                                  init_schedule.blocks[b-1].departure_time;

      double dts = dt[randomIndex]*largest_step_size; 
      // apply move to schedule blocks

      if ((dts < 0)|| (std::isnan(dts)) || (std::isinf(dts))) {

        list_of_schedules.push_back(schedule_sol);
        deltaVs_of_neighbourhood.push_back(std::numeric_limits<double>::infinity());
        
        continue;
      }
      move_wrapper(schedule_sol.blocks, b, dts, move_current,simfile);


      double DeltaVMinimaopt = 0.0;

      // Quick feasibility check with new constraint logic
      if (!quick_feasibility_check(schedule_sol, b)){

        list_of_schedules.push_back(schedule_sol);
        deltaVs_of_neighbourhood.push_back(std::numeric_limits<double>::infinity());
        continue;
      
      }

      // Process arrival moves (mapped to "add"/"sub")
      if ((move_current == "add" )|| (move_current == "sub"))
      {
          find_optimal_trajectory_no_iter(
              schedule_sol.blocks[b - 1].satname,
              schedule_sol.blocks[b].satname,
              schedule_sol.blocks[b - 1].departure_time,
              schedule_sol.blocks[b].arrival_time, simfile, DeltaVMinimaopt);

          schedule_sol.blocks[b].deltaV_arrival = DeltaVMinimaopt;
        
      }

      // Process swap moves with first/last block protection
      if (((move_current == "swap") || (move_current == "swap_inv")) && 
          (b+1 < init_schedule.blocks.size()-1))
      {
        // if move did nothing
        if (schedule_sol.blocks[b].satname == init_schedule.blocks[b].satname){
          
          list_of_schedules.push_back(schedule_sol);
          deltaVs_of_neighbourhood.push_back(std::numeric_limits<double>::infinity());
          continue;
        }

        // Swap moves already handle deltaV recalculation in swap_slots functions
        // No additional deltaV calculation needed for performance
      }
      list_of_schedules.push_back(schedule_sol);

      // closes for loop over the whole schedule

      double totalDeltaV_of_sol = 0.0;

      for (int elem = 0; elem < schedule_sol.blocks.size(); elem++)
      {
        totalDeltaV_of_sol += std::abs(schedule_sol.blocks[elem].deltaV_arrival);
      }

      deltaVs_of_neighbourhood.push_back(totalDeltaV_of_sol);

      if (totalDeltaV_of_sol == 0)
      {
        std::cout << "0 delta V";
      }

    } // closes iteration over the full schedule

    if (deltaVs_of_neighbourhood.empty()) {
      move_choice += 1;
      if (move_choice -1 >= (move_methods.size())) {
        move_choice -= move_choice;
        move_choice += 1;
      }
      continue;
    }
    neighbourhood_minima = find_minima_val(deltaVs_of_neighbourhood);
    // std::cout<< "neighbourhood_minima : "<< neighbourhood_minima <<
    // std::endl; std::cout<< "neighbourhood_maxima : "<<
    // deltaVs_of_neighbourhood.max() << std::endl; std::cout<< "size : " <<
    // deltaVs_of_neighbourhood.size()<< std::endl;

    // if no improvement is found, then stop
    if (neighbourhood_minima >= deltaVminima_so_far)
    {
      move_choice += 1;
      if (move_choice -1 >= (move_methods.size()))
      {
        move_choice -= move_choice;
        move_choice += 1;
      }


      int index_minima = find_minima_index(deltaVs_of_neighbourhood);

    }

    else
    {
      int index_minima = find_minima_index(deltaVs_of_neighbourhood);

      deltaVminima_so_far = std::abs(neighbourhood_minima);
      init_schedule = list_of_schedules[index_minima];
      init_deltaV = deltaVminima_so_far;
      // view_schedule(list_of_schedules[index_minima]);
    }

  } // closes while
}

void run_vn_search(DataFrame simfile, std::vector<double> move_size,
                   std::vector<std::string> moves_to_consider,
                   std::vector<std::string> sat_names_in_schedule,
                   std::vector<double> t_depart, std::vector<double> t_arrive,
                   double &deltaV_of_schedule, double service_time, int max_iter)
{

  // construct initial schedule using provided departure and arrival times
  schedule_struct init_schedule = create_schedule_lambert_only(deltaV_of_schedule, t_arrive, t_depart,
                                                               sat_names_in_schedule, simfile, service_time);

  std::cout << "\n";

  std::cout << "Initial Schedule" << std::endl;

  view_schedule(init_schedule);

  std::cout << "\n";

  std::cout << "Total Delta V: " << deltaV_of_schedule;

  // Now the optimal schedule is calculated using local search
  vn_search(deltaV_of_schedule, init_schedule, move_size, simfile, service_time, moves_to_consider, max_iter);

  std::cout << "\n";

  std::cout << "Result of Local Search" << std::endl;

  view_schedule(init_schedule);

  std::cout << "\n";

  std::cout << "Total Delta V: " << deltaV_of_schedule << "\n";
}

void run_vn_search_fixed_tarrive(DataFrame simfile, std::vector<double> move_size,
                                 std::vector<std::string> moves_to_consider,
                                 std::vector<std::string> sat_names_in_schedule,
                                 std::vector<double> t_depart, std::vector<double> t_arrive,
                                 double &deltaV_of_schedule, double service_time, int max_iter)
{

  // construct initial schedule using provided departure and arrival times
  schedule_struct init_schedule = create_schedule_lambert_only(deltaV_of_schedule, t_arrive, t_depart, sat_names_in_schedule, simfile, service_time);

  // divide the schedule into chunks of 2 and run vns on each block
  double pair_deltaV = 0;

  for (int b = 1; b < init_schedule.blocks.size(); b++)
  {

    schedule_struct pair_to_pass;
    pair_to_pass.blocks.push_back(init_schedule.blocks[b - 1]);
    pair_to_pass.blocks.push_back(init_schedule.blocks[b]);

    pair_deltaV = init_schedule.blocks[b - 1].deltaV_arrival + init_schedule.blocks[b].deltaV_arrival;

    // since we will have n-1 transfers to consider for n blocks we give each transfer max_iter/(size-1)
    vn_search(pair_deltaV, pair_to_pass, move_size, simfile, service_time, moves_to_consider, (max_iter / ((int)init_schedule.blocks.size() - 1)));

    init_schedule.blocks[b - 1] = pair_to_pass.blocks[0];
    init_schedule.blocks[b] = pair_to_pass.blocks[1];
  }

  double deltav_of_full_schedule = 0;

  for (task_block bl : init_schedule.blocks)
  {

    deltav_of_full_schedule += bl.deltaV_arrival;
  }

  deltaV_of_schedule = deltav_of_full_schedule;

  view_schedule(init_schedule);

  std::cout << "Total Delta V: " << deltaV_of_schedule << "\n";
}
