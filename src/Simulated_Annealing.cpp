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
#include "Local_search.hpp"
#include <random>
/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/

double prob_calculation(double old_state_energy, double new_state_energy,double T){

      return std::pow(2,-(new_state_energy - old_state_energy)/T); 

}



schedule_struct simulated_annealing_lambert(double &init_deltaV, schedule_struct init_schedule, std::vector<double> dt_move,
                                                        DataFrame simfile, double service_time,std::vector<std::string> move_methods, double cooling_param, int maxiter, double Temp){

  double deltaVminima_so_far = init_deltaV;

  schedule_struct optimal_schedule; 

  int iternum = 0 ;  

  std::vector<double> arrival_constraints; 
  std::vector<double> departure_constraints; 
  for(int b=0 ; b < init_schedule.blocks.size(); b++){

    arrival_constraints.push_back(init_schedule.blocks[b].arrival_time);
    departure_constraints.push_back(init_schedule.blocks[b].departure_time);
  }

  std::random_device rd;  // Used to seed the generator
  std::mt19937 gen(rd()); // Mersenne Twister generator
  std::uniform_real_distribution<> dis(0.0, 1.0);


  while (iternum < maxiter){
    
    iternum +=1; 
    Temp = Temp*cooling_param;  
    std::cout << "Temp:" << Temp<< std::endl;
    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;

    // compute neighbourhood size safely to avoid integer overflow or
    // attempts to allocate an excessively large vector which throws
    // std::bad_array_new_length
    size_t neighbourhood_size = init_schedule.blocks.size();
    neighbourhood_size *= move_methods.size();
    neighbourhood_size *= dt_move.size();

    if (neighbourhood_size == 0) {
      // nothing to explore
      break;
    }

    const size_t MAX_NEIGHBOURHOOD = 10000000; // sanity cap
    if (neighbourhood_size > MAX_NEIGHBOURHOOD) {
      std::cerr << "Neighbourhood size too large: " << neighbourhood_size << ". Aborting search to avoid OOM." << std::endl;
      break;
    }

    arma::vec deltaVs_of_neighbourhood(static_cast<arma::uword>(neighbourhood_size));

    double neighbourhood_minima;

    int solnum_neighbourhood = 0;


    //loop over all moves 
    for (int m = 0; m < move_methods.size(); m++){
      for (int d = 0; d < dt_move.size(); d++){
        for (int b = 0; b < init_schedule.blocks.size(); b++){

        schedule_struct schedule_sol = init_schedule;
        
        // apply move to schedule blocks
        move_wrapper(schedule_sol.blocks, b, dt_move[d], move_methods[m],
                       simfile);
        
        
        double DeltaVMinimaopt;

        if((move_methods[m] == "add departure") || (move_methods[m]== "sub departure")){
              
            bool b_not_last_elem = b < (schedule_sol.blocks.size() - 1);

            if (b_not_last_elem){
              
              bool condition_add = ((move_methods[m] == "add departure") && (schedule_sol.blocks[b].departure_time < schedule_sol.blocks[b+1].arrival_time));
              bool condition_sub = ((move_methods[m] == "sub departure") && (departure_constraints[b] < schedule_sol.blocks[b].departure_time));
              

              if (condition_add || condition_sub){
              
                //finding lowest energy lambert transfer for given tof
                find_optimal_trajectory_no_iter(
                  schedule_sol.blocks[b].satname,
                  schedule_sol.blocks[b + 1].satname,
                  schedule_sol.blocks[b].departure_time,
                  schedule_sol.blocks[b + 1].arrival_time, simfile, DeltaVMinimaopt);

              
                schedule_sol.blocks[b+1].deltaV_arrival = DeltaVMinimaopt;

              
              }  
            }
          }


        if((move_methods[m] == "add arrival") || (move_methods[m]== "sub arrival")){
        
            if(b > 0){

              bool condition_add = ((move_methods[m] == "add arrival") && (schedule_sol.blocks[b].arrival_time < arrival_constraints[b]));
              bool condition_sub = ((move_methods[m] == "sub arrival") && (schedule_sol.blocks[b-1].departure_time < schedule_sol.blocks[b].arrival_time));


              if (condition_add || condition_sub) {


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
        
        //closes for loop over the whole schedule

        double totalDeltaV_of_sol = 0.0;

        for (int elem = 0; elem < schedule_sol.blocks.size(); elem++) {

          totalDeltaV_of_sol +=
              std::abs(schedule_sol.blocks[elem].deltaV_arrival);
        }

        deltaVs_of_neighbourhood(solnum_neighbourhood) = (totalDeltaV_of_sol);

        if (totalDeltaV_of_sol == 0) {
          std::cout << "0 delta V";
        }

        solnum_neighbourhood += 1;

        } // closes iteration over the full schedule 
      }//closes iteration over all move sizes
    } // closes iteration over all the moves 

    
    neighbourhood_minima = deltaVs_of_neighbourhood.min();
    // std::cout<< "neighbourhood_minima : "<< neighbourhood_minima <<
    // std::endl; std::cout<< "neighbourhood_maxima : "<<
    // deltaVs_of_neighbourhood.max() << std::endl; std::cout<< "size : " <<
    // deltaVs_of_neighbourhood.size()<< std::endl;

    // if no improvement is found, 
    if (neighbourhood_minima > deltaVminima_so_far) {
     
      double randomValue = dis(gen); 
      
      std::cout<<"prob calc:"<< prob_calculation(deltaVminima_so_far,neighbourhood_minima,Temp) <<std::endl;
      std::cout<<"random val:"<< randomValue<<std::endl;
      std::cout<<"Delta V minima diff:"<< deltaVminima_so_far - neighbourhood_minima <<"\n"; 

      if (prob_calculation(deltaVminima_so_far,neighbourhood_minima,Temp) > randomValue){

        init_deltaV = neighbourhood_minima;
        arma::uword index_minima_uword = arma::index_min(deltaVs_of_neighbourhood);
        int index_minima = static_cast<int>(index_minima_uword);
        deltaVminima_so_far = std::abs(neighbourhood_minima);
        init_schedule = list_of_schedules[index_minima];
      
      }
      
      // view_schedule(list_of_schedules[index_minima]);
    }


    else if(neighbourhood_minima <= deltaVminima_so_far){

      init_deltaV = neighbourhood_minima;
      arma::uword index_minima_uword = arma::index_min(deltaVs_of_neighbourhood);
      int index_minima = static_cast<int>(index_minima_uword);
      deltaVminima_so_far = std::abs(neighbourhood_minima);
      init_schedule = list_of_schedules[index_minima];

    }


  }//closes while
  
  optimal_schedule = init_schedule;
  return optimal_schedule;


}




void run_simulated_annealing( DataFrame simfile, std::vector<double> move_size,  
                      std::vector<std::string> moves_to_consider,
                      std::vector<std::string> sat_names_in_schedule ,
                      std::vector<double> t_depart, std::vector<double> t_arrive, 
                      double &deltaV_of_schedule, double service_time, double cooling_param, int maxiter, double T_init){


    //construct initial schedule using provided departure and arrival times 
    schedule_struct init_schedule = create_schedule_lambert_only(deltaV_of_schedule, t_arrive, t_depart, 
                                                                  sat_names_in_schedule,simfile, service_time);

    std::cout<< "\n"; 

    std::cout << "Initial Schedule"<< std::endl; 

    view_schedule(init_schedule);
    
    std::cout<< "\n"; 

    std::cout << "Total Delta V: "<< deltaV_of_schedule;

    // Now the optimal schedule is calculated using local search
    schedule_struct findopt_schedule = simulated_annealing_lambert( deltaV_of_schedule, init_schedule, move_size, 
                                                                    simfile, service_time, moves_to_consider, cooling_param, maxiter, T_init);
    
    std::cout<<"\n";

    std::cout << "Result of Local Search"<<std::endl; 

    view_schedule(findopt_schedule);

    std::cout<<"\n";

    std::cout<<"Total Delta V: "<< deltaV_of_schedule<<"\n";


    


}

































