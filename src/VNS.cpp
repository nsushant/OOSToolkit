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

double find_minima_val(std::vector<double> v)
{

  auto itmin = std::min_element(v.begin(), v.end());
  double minima = *itmin; // smallest element from dereferenced iterator

  return minima;
}

int find_minima_index(std::vector<double> v)
{

  auto itmin = std::min_element(v.begin(), v.end());
  int index_minima = std::distance(v.begin(), itmin);

  return index_minima;
}

void vn_search(double &init_deltaV, schedule_struct& init_schedule, std::vector<double> dt_move,
                          DataFrame simfile, double service_time, std::vector<std::string> move_methods, int max_iter){

  double temp = max_iter/2;


  std::map<int, std::pair<std::string, double>> possible_moves;
  int k = 0;
  int iterval = 0;

  for (int mov = 0; mov < move_methods.size(); mov++)
  {
    for (int dts = 0; dts < dt_move.size(); dts++)
    {

      std::pair<std::string, double> keyp = {move_methods[mov], dt_move[dts]};

      possible_moves[k] = keyp;
      k += 1;
    }
  }

  double deltaVminima_so_far = init_deltaV;
  schedule_struct optimal_schedule;

  std::vector<double> arrival_constraints;
  std::vector<double> departure_constraints;

  for (int b = 0; b < init_schedule.blocks.size(); b++)
  {
    arrival_constraints.push_back(init_schedule.blocks[b].arrival_time);
    departure_constraints.push_back(init_schedule.blocks[b].departure_time);
  }

  int move_choice = 0;



  std::random_device rd;                      // non-deterministic seed
  std::mt19937 gen(rd());                     // Mersenne Twister RNG
                                              //
  std::uniform_real_distribution<> dist(0.0, 1.0); // range [1, 100]

  double decay = 0.9; 

  while (iterval < max_iter)
  { temp*=0.6;
    iterval +=1;
    std::pair<std::string, double> move_current = possible_moves[move_choice];
    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;
    std::vector<double> deltaVs_of_neighbourhood;

    double neighbourhood_minima;

    // loop over all moves
    for (int b = 0; b < init_schedule.blocks.size(); b++)
    {

      schedule_struct schedule_sol = init_schedule;

      // apply move to schedule blocks
      move_wrapper(schedule_sol.blocks, b, move_current.second, move_current.first,
                   simfile);

      double DeltaVMinimaopt;

      if ((move_current.first == "add departure") || (move_current.first == "sub departure"))
      {

        bool b_not_last_elem = b < (schedule_sol.blocks.size() - 1);

        if (b_not_last_elem)
        {

          bool condition_add = ((move_current.first == "add departure") && (schedule_sol.blocks[b].departure_time < schedule_sol.blocks[b + 1].arrival_time));
          bool condition_sub = ((move_current.first == "sub departure") && (departure_constraints[b] <= schedule_sol.blocks[b].departure_time));

          if (condition_add || condition_sub)
          {

            // finding lowest energy lambert transfer for given tof
            find_optimal_trajectory_no_iter(
                schedule_sol.blocks[b].satname,
                schedule_sol.blocks[b + 1].satname,
                schedule_sol.blocks[b].departure_time,
                schedule_sol.blocks[b + 1].arrival_time, simfile, DeltaVMinimaopt);

            schedule_sol.blocks[b + 1].deltaV_arrival = DeltaVMinimaopt;
          }
        }
      }

      if ((move_current.first == "add arrival") || (move_current.first == "sub arrival"))
      {

        if (b > 0)
        {

          bool condition_add = ((move_current.first == "add arrival") && (schedule_sol.blocks[b].arrival_time <= arrival_constraints[b]));
          bool condition_sub = ((move_current.first == "sub arrival") && (schedule_sol.blocks[b - 1].departure_time < schedule_sol.blocks[b].arrival_time));

          if (condition_add || condition_sub)
          {

            find_optimal_trajectory_no_iter(
                schedule_sol.blocks[b - 1].satname,
                schedule_sol.blocks[b].satname,
                schedule_sol.blocks[b - 1].departure_time,
                schedule_sol.blocks[b].arrival_time, simfile, DeltaVMinimaopt);

            schedule_sol.blocks[b].deltaV_arrival = DeltaVMinimaopt;
          }
        }
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

    neighbourhood_minima = find_minima_val(deltaVs_of_neighbourhood);
    // std::cout<< "neighbourhood_minima : "<< neighbourhood_minima <<
    // std::endl; std::cout<< "neighbourhood_maxima : "<<
    // deltaVs_of_neighbourhood.max() << std::endl; std::cout<< "size : " <<
    // deltaVs_of_neighbourhood.size()<< std::endl;

    // if no improvement is found, then stop
    if (neighbourhood_minima >= deltaVminima_so_far)
    {
        
        
        move_choice += 4; 
        /*
        //std::cout <<"step = "<<random_step << std::endl;

        // we need to wrap around the max array len 
        

        if(temp <= (random_step/(double)possible_moves.size()) ){

        move_choice += random_step; 
        move_choice = move_choice % (possible_moves.size()-1);
        
        }*/
        
        double move_random = dist(gen); 

        decay *=0.70;

        if (move_random < decay){
        
        int index_minima = find_minima_index(deltaVs_of_neighbourhood);

        deltaVminima_so_far = std::abs(neighbourhood_minima);
        init_schedule = list_of_schedules[index_minima];

        //move_choice +=1;
        
        } 
           

        
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
