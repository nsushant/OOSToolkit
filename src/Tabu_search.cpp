#include "Local_search.hpp"
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"
#include <armadillo>
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <map>
#include <cassert>
#include "Local_search.hpp"

/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/

// moves
// use local search to find the optimal_schedule

schedule_struct Tabu_search_opt_schedule_lambert_only(double &init_deltaV, schedule_struct &init_schedule, std::vector<double> dt_move,
                                                      DataFrame simfile, double service_time, std::vector<std::string> move_methods)
{

  int saturation_limit = 10;
  int maxiter = 50;
  double deltaVminima_so_far = init_deltaV;

  std::map<std::pair<std::string, double>, int> possible_moves;

  for (int mov = 0; mov < move_methods.size(); mov++)
  {

    for (int dts = 0; dts < dt_move.size(); dts++)
    {

      std::pair<std::string, double> keyp = {move_methods[mov], dt_move[dts]};

      possible_moves[keyp] = 0;
    }
  }

  schedule_struct optimal_schedule;

  std::vector<double> arrival_constraints;
  std::vector<double> departure_constraints;

  for (int b = 0; b < init_schedule.blocks.size(); b++)
  {

    arrival_constraints.push_back(init_schedule.blocks[b].arrival_time);
    departure_constraints.push_back(init_schedule.blocks[b].departure_time);
  }

  std::map<std::pair<std::string, double>, int> tabu_moves;

  int iternum = 0;

  while (iternum < maxiter)
  {
    iternum += 1;
    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;
    std::vector<std::string> move_chosen;
    std::vector<double> move_size_chosen;

    // compute neighbourhood size safely to avoid integer overflow or
    // attempts to allocate an excessively large vector which throws
    // std::bad_array_new_length

    std::vector<double> deltaVs_of_neighbourhood;
    bool stop_all = false;

    double neighbourhood_minima;

    int solnum_neighbourhood = 0;
    std::cout << "startloop" << std::endl;
    // loop over all moves
    for (int m = 0; m < move_methods.size() && !stop_all; m++)
    {
      for (int d = 0; d < dt_move.size() && !stop_all; d++)
      {
        for (int b = 0; b < init_schedule.blocks.size(); b++)
        {

          std::pair<std::string, double> keypair = {move_methods[m], dt_move[d]};
          bool move_is_tabu = ((tabu_moves.find(keypair) != tabu_moves.end()));

          schedule_struct schedule_sol = init_schedule;

          // apply move to schedule blocks
          move_wrapper(schedule_sol.blocks, b, dt_move[d], move_methods[m],
                       simfile);

          double DeltaVMinimaopt;

          if ((move_methods[m] == "add departure") || (move_methods[m] == "sub departure"))
          {

            bool b_not_last_elem = b < (schedule_sol.blocks.size() - 1);

            if (b_not_last_elem)
            {

              bool condition_add = ((move_methods[m] == "add departure") && (schedule_sol.blocks[b].departure_time <= schedule_sol.blocks[b + 1].arrival_time));
              bool condition_sub = ((move_methods[m] == "sub departure") && (departure_constraints[b] < schedule_sol.blocks[b].departure_time));

              if (condition_add || condition_sub)
              {

                // finding lowest energy lambert transfer for given tof
                find_optimal_trajectory_no_iter(
                    schedule_sol.blocks[b].satname,
                    schedule_sol.blocks[b + 1].satname,
                    schedule_sol.blocks[b].departure_time,
                    schedule_sol.blocks[b + 1].arrival_time, simfile, DeltaVMinimaopt);

                schedule_sol.blocks[b + 1].deltaV_arrival = DeltaVMinimaopt;
                // move_chosen.push_back(move_methods[m]);
                // move_size_chosen.push_back(dt_move[d]);
              }
            }
          }

          if ((move_methods[m] == "add arrival") || (move_methods[m] == "sub arrival"))
          {

            if (b > 0)
            {

              bool condition_add = ((move_methods[m] == "add arrival") && (schedule_sol.blocks[b].arrival_time <= arrival_constraints[b]));
              bool condition_sub = ((move_methods[m] == "sub arrival") && (schedule_sol.blocks[b - 1].departure_time < schedule_sol.blocks[b].arrival_time));

              if (condition_add || condition_sub)
              {

                find_optimal_trajectory_no_iter(
                    schedule_sol.blocks[b - 1].satname,
                    schedule_sol.blocks[b].satname,
                    schedule_sol.blocks[b - 1].departure_time,
                    schedule_sol.blocks[b].arrival_time, simfile, DeltaVMinimaopt);

                schedule_sol.blocks[b].deltaV_arrival = DeltaVMinimaopt;

                // move_chosen.push_back(move_methods[m]);
                // move_size_chosen.push_back(dt_move[d]);
              }
            }
          }

          double totalDeltaV_of_sol = 0.0;

          for (int elem = 0; elem < schedule_sol.blocks.size(); elem++)
          {

            totalDeltaV_of_sol +=
                std::abs(schedule_sol.blocks[elem].deltaV_arrival);
          }

          if (move_is_tabu)
          {
            std::cout << "Move is tabu, skipping " << std::endl;
            tabu_moves[keypair] -= 1;
            std::cout << "value after increment" << tabu_moves[keypair] << "\n";

            if (tabu_moves[keypair] == 0)
            {
              tabu_moves.erase(keypair);
              std::cout << "erased keypair from tabu list" << std::endl;
            }

            if (deltaVminima_so_far > totalDeltaV_of_sol)
            {
              std::cout << "skipping" << std::endl;
              continue;
            }
          }

          move_chosen.push_back(move_methods[m]);
          move_size_chosen.push_back(dt_move[d]);
          list_of_schedules.push_back(schedule_sol);

          // closes for loop over the whole schedule

          deltaVs_of_neighbourhood.push_back(totalDeltaV_of_sol);

          if (totalDeltaV_of_sol == 0)
          {
            std::cout << "0 delta V";
          }

          solnum_neighbourhood += 1;

        } // closes iteration over the full schedule
      } // closes iteration over all move sizes
    } // closes iteration over all the moves

    auto itminDeltaV = std::min_element(deltaVs_of_neighbourhood.begin(), deltaVs_of_neighbourhood.end());
    double deltaVs_of_neighbourhood_minima = *itminDeltaV; // smallest element from dereferenced iterator

    neighbourhood_minima = deltaVs_of_neighbourhood_minima;

    // if no improvement is found, then stop
    if (neighbourhood_minima >= deltaVminima_so_far)
    {
      if (neighbourhood_minima == deltaVminima_so_far)
      {

        if (saturation_limit == 0)
        {
          std::cout << "reached saturation limit" << std::endl;
          break;
        }

        saturation_limit -= 1;
      }

      std::size_t num_tabu = tabu_moves.size();
      std::size_t num_total_m = possible_moves.size();

      if (num_total_m == num_tabu)
      {

        std::cout << "All moves Tabu, stopping search" << num_total_m << "," << num_tabu << "\n";
        break;
      }

      std::cout << "Stuck in minima, assigning tabu" << std::endl;

      int index_minima = std::distance(deltaVs_of_neighbourhood.begin(), itminDeltaV);
      deltaVminima_so_far = std::abs(neighbourhood_minima);

      int midx_minima = index_minima;
      int didx_minima = index_minima;

      init_schedule = list_of_schedules[index_minima];

      assert(deltaVs_of_neighbourhood.size() == move_chosen.size());
      std::pair<std::string, double> tabu_keypair = {move_chosen[midx_minima], move_size_chosen[didx_minima]};

      std::cout << "created keypair" << std::endl;

      tabu_moves[tabu_keypair] = 30;

      std::cout << "created tabu move" << std::endl;
    }

    // else if improvement is found, adopt new solution
    else
    {
      std::cout << "Improvement found" << std::endl;
      int index_minima = std::distance(deltaVs_of_neighbourhood.begin(), itminDeltaV);

      deltaVminima_so_far = std::abs(neighbourhood_minima);
      init_schedule = list_of_schedules[index_minima];
      init_deltaV = deltaVminima_so_far;
      // view_schedule(list_of_schedules[index_minima]);
    }

  } // closes while

  return init_schedule;
}

void run_tabu_search(DataFrame simfile, std::vector<double> move_size,
                     std::vector<std::string> moves_to_consider,
                     std::vector<std::string> sat_names_in_schedule,
                     std::vector<double> t_depart, std::vector<double> t_arrive,
                     double &deltaV_of_schedule, double service_time)
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
  schedule_struct findopt_schedule = Tabu_search_opt_schedule_lambert_only(deltaV_of_schedule, init_schedule, move_size,
                                                                           simfile, service_time, moves_to_consider);

  std::cout << "\n";

  std::cout << "Result of Local Search" << std::endl;

  view_schedule(findopt_schedule);

  std::cout << "\n";

  std::cout << "Total Delta V: " << deltaV_of_schedule << "\n";
}
