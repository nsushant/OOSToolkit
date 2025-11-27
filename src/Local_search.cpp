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
/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local
search.

*/

// moves
void move_dt(task_block &tb, double dt) {

  // requires both lambert tranfers to and from the block to be 
  // recomputed

  tb.arrival_time -= dt;
  tb.departure_time -= dt;

}



void move_dt2(std::vector<task_block> &blocks, int b_index, double dt) {

  if (b_index == 0) {

    blocks[b_index].departure_time += dt; /// 2 - ((int)dt / 2 % 100);

  }

  else {

    if (blocks[b_index].arrival_time > blocks[b_index - 1].departure_time) {

      blocks[b_index].arrival_time -= dt;
      // blocks[b_index].arrival_time += dt;

      if (blocks[b_index].departure_time == 0) {

        blocks[b_index].departure_time += 0;

      }

      else {

        blocks[b_index].departure_time += dt / 2 - ((int)dt / 2 % 100);
      }
    }
  }
}
void move_add_arrival(std::vector<task_block> &blocks, int b_index, double dt) {

    //except for the first block
    if (b_index > 0 ) {
        
    // ensure that the addition in arrival time does not make arrival occur after the departure time
    if ((blocks[b_index].arrival_time + blocks[b_index].service_duration + dt) < blocks[b_index].departure_time) {

        blocks[b_index].arrival_time += dt;

    }


    // since the departure time of the last block is 0 

    else if(b_index == (blocks.size()-1)){blocks[b_index].arrival_time += dt;}

    }


}
void move_add_departure(std::vector<task_block> &blocks, int b_index,
                        double dt) {

  //for all but the second to last block                          
  if ( b_index != (blocks.size() - 1) ) {

    // make sure the shuttle departs before the next expected arrival
    if (blocks[b_index].departure_time + dt <
        blocks[b_index + 1].arrival_time) {

      blocks[b_index].departure_time += dt;
    }


  }


}

void move_sub_departure(std::vector<task_block> &blocks, int b_index,
                        double dt) {
  
  // for all but the last block
  if (b_index != (blocks.size() - 1)) {

    if (blocks[b_index].departure_time - dt > (blocks[b_index].arrival_time + blocks[b_index].service_duration)) {

      blocks[b_index].departure_time -= dt;
    }
  }



}


void move_sub_arrival(std::vector<task_block> &blocks, int b_index, double dt){

  // for all but the first block
  if (b_index > 0 ) {

    if (blocks[b_index].arrival_time - dt > blocks[b_index - 1].departure_time) {

      blocks[b_index].arrival_time -= dt;
    }

  }


}


void swap_slots(std::vector<task_block> &blocks, int b_index, double dt,
                DataFrame simfile) {

  double deadline = blocks[b_index].demand_deadline;

  // get valid blocks
  std::vector<int> candidate_blocks;

  if (b_index % 2 == 0) {

    task_block copyblock = blocks[b_index + 2];

    if (copyblock.departure_time < deadline) {

      double service_minimum =
          blocks[b_index].arrival_time + blocks[b_index].service_duration;

      blocks[b_index + 2].satname = blocks[b_index].satname;
      blocks[b_index + 2].service_duration = blocks[b_index].service_duration;
      blocks[b_index + 2].demand_deadline = blocks[b_index].demand_deadline;

      blocks[b_index].demand_deadline = copyblock.demand_deadline;
      blocks[b_index].service_duration = copyblock.service_duration;
      blocks[b_index].satname = copyblock.satname;

      bool timeslot_not_long_enough = (blocks[b_index + 2].departure_time <
                                       (blocks[b_index + 2].arrival_time +
                                        blocks[b_index + 2].service_duration));

      if (timeslot_not_long_enough) {

        double time_adjustment_needed = (blocks[b_index + 2].arrival_time +
                                         blocks[b_index + 2].service_duration) -
                                        blocks[b_index + 2].departure_time;

        blocks[b_index + 2].arrival_time -= time_adjustment_needed;

        blocks[b_index + 1].departure_time -= time_adjustment_needed;
        blocks[b_index + 1].arrival_time -= time_adjustment_needed;

        blocks[b_index].departure_time -= time_adjustment_needed;

        std::cout << "proposed switch was invalid, departure time < "
                     "service_time + arrival_time, adjustment made";
      }

      bool timeslot_not_long_enough_2 =
          (blocks[b_index].departure_time <
           (blocks[b_index].arrival_time + blocks[b_index].service_duration));

      if (timeslot_not_long_enough_2) {

        double time_adjustment_needed =
            (blocks[b_index].arrival_time + blocks[b_index].service_duration) -
            blocks[b_index].departure_time;

        blocks[b_index].departure_time += time_adjustment_needed;

        blocks[b_index + 1].departure_time += time_adjustment_needed;
        blocks[b_index + 1].arrival_time += time_adjustment_needed;

        blocks[b_index + 2].arrival_time += time_adjustment_needed;
      }

      // recalculate delta Vs
      find_optimal_trajectory_no_iter(
          blocks[b_index - 1].satname, blocks[b_index].satname,
          blocks[b_index - 1].departure_time, blocks[b_index].arrival_time,
          simfile, blocks[b_index].deltaV_arrival);
      find_optimal_trajectory_no_iter(
          blocks[b_index].satname, blocks[b_index + 1].satname,
          blocks[b_index].departure_time, blocks[b_index + 1].arrival_time,
          simfile, blocks[b_index + 1].deltaV_arrival);
      find_optimal_trajectory_no_iter(
          blocks[b_index + 1].satname, blocks[b_index + 2].satname,
          blocks[b_index + 1].departure_time, blocks[b_index + 2].arrival_time,
          simfile, blocks[b_index + 2].deltaV_arrival);
      find_optimal_trajectory_no_iter(
          blocks[b_index + 2].satname, blocks[b_index + 3].satname,
          blocks[b_index + 2].departure_time, blocks[b_index + 3].arrival_time,
          simfile, blocks[b_index + 3].deltaV_arrival);
    }
  }
}

void move_wrapper(std::vector<task_block> &blocks, int b_index, double dt,
                  std::string method, DataFrame simfile) {

  if (method == "add arrival") {
    move_add_arrival(blocks, b_index, dt);
  }
  if (method == "add departure") {
    move_add_departure(blocks, b_index, dt);
  }
  if (method == "sub arrival") {
    move_sub_arrival(blocks, b_index, dt);
  }

  if (method == "sub departure") {
    move_sub_departure(blocks, b_index, dt);
  }
  
  if (method == "move_dt2") {
    move_dt2(blocks, b_index, dt);
  }
  
  if (method == "swap_slots") {
    swap_slots(blocks, b_index, dt, simfile);
  }
}

// allows you to print out the schedule struct in a readable format
void view_schedule(schedule_struct schedule_to_print) {

  // Define column widths
  const int w1 = 15; // Task
  const int w2 = 15; // Arrival
  const int w3 = 15; // Departure
  const int w4 = 13; // Fuel cost

  // Print header line
  std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";
  std::cout << "| " << std::left << std::setw(w1 - 1) << "Task"
            << "| " << std::setw(w2 - 1) << "Arrival Time"
            << "| " << std::setw(w3 - 1) << "Departure Time"
            << "| " << std::setw(w4 - 1) << "Delta V"
            << "|\n";
  std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";

  // Print each task
  for (const auto &t : schedule_to_print.blocks) {
    std::cout << "| " << std::left << std::setw(w1 - 1) << t.satname << "| "
              << std::setw(w2 - 1) << t.arrival_time << "| "
              << std::setw(w3 - 1) << t.departure_time << "| " << std::right
              << std::setw(w4 - 3) << std::fixed << std::setprecision(2)
              << t.deltaV_arrival << "  |\n";
  }

  std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";
}


void init_dep_arrival_times_random(std::vector<double> &departure_times, std::vector<double> &arrival_times,double period, double service_time, int num_sats){
        
        // if we have n sats to visit, we will have 2n +1 task blocks

        int num_blocks = 2*num_sats + 1;

        departure_times.push_back(0.0);
        arrival_times.push_back(0.0);
        
        //seed random number generator
        std::random_device rdom;  
        std::mt19937 gen(rdom());
        std::uniform_int_distribution<> distrib(1,3);

        for(int i = 1 ; i < num_blocks ; i++){
         
          int orbmul = distrib(gen);
          double orbt = period*orbmul;
          
          arrival_times.push_back(departure_times[i-1] + orbt); 
          departure_times.push_back(arrival_times[i] + service_time);
          
        }

}


void init_dep_arrival_times_strict_timespan(std::vector<double> &departure_times, std::vector<double> &arrival_times, double final_time, double service_time, int num_sats){
      
      int num_blocks = 2*num_sats + 1;

      departure_times.push_back(0.0); 
      arrival_times.push_back(0.0); 
      
      int dt = ((int)final_time/num_blocks) - ((int)final_time/num_blocks % 100); 

      double t_tot = 0.0;

      for (int i = 1; i < num_blocks; i++){
      
        t_tot+=dt; 
        arrival_times.push_back(t_tot);
        departure_times.push_back(t_tot + service_time); 
        

        if (arrival_times[i] > departure_times[i] + dt){
        
          std::cout<<"May contain overlapping blocks"<<std::endl; 
        
        }

      }

}


std::vector<std::string> init_satname_array(std::string service_satname, std::vector<std::string> client_satnames, bool allow_revisits, double num_sats){

    int num_blocks = 2*num_sats + 1; 

    std::vector<std::string> satnames; 

    satnames.push_back(service_satname); 

    int d = 0;

    for (int i = 1; i < num_blocks; i++){

        if (i % 2 == 0){
          
          satnames.push_back(service_satname); 
          
        }
          

        else{

          if (allow_revisits == true){
            
            if(d==client_satnames.size()){

              d-=d; 

            }
            
            satnames.push_back(client_satnames[d]);
            d+=1;

          }

          else{

            satnames.push_back(client_satnames[d]);
            d+=1;
            
          }

        }


    }
    
    return satnames; 
} 


// initialize a schedule
schedule_struct create_schedule(double &deltaV_of_schedule_init,
                                std::vector<double> arrival_times,
                                std::vector<double> departure_times,
                                std::vector<std::string> satnames,
                                DataFrame simfile) {

  // initial block for the service depot

  schedule_struct schedule;

  task_block initblock;

  initblock.arrival_time = arrival_times[0];

  initblock.departure_time = departure_times[0];

  initblock.satname = satnames[0];

  initblock.deltaV_arrival = 0.0;

  schedule.blocks.push_back(initblock);

  deltaV_of_schedule_init = 0.0;

  for (int i = 1; i < satnames.size(); i++) {

    task_block block;

    block.arrival_time = arrival_times[i];

    block.departure_time = departure_times[i];

    block.satname = satnames[i];

    double deltaV_transfer = compact_optimal_calc(
        satnames[i - 1], block.satname, departure_times[i - 1],
        block.arrival_time, simfile);

    block.deltaV_arrival = deltaV_transfer;

    deltaV_of_schedule_init += deltaV_transfer;

    schedule.blocks.push_back(block);
  }

  return schedule;
}

schedule_struct create_schedule_lambert_only(
    double &deltaV_of_schedule_init, std::vector<double> arrival_times,
    std::vector<double> departure_times, std::vector<std::string> satnames,
    DataFrame simfile, double service_time, std::vector<double> deadlines){

  // initial block for the service depot

  schedule_struct schedule;

  task_block initblock;

  initblock.arrival_time = arrival_times[0];

  initblock.departure_time = departure_times[0];

  initblock.satname = satnames[0];

  initblock.deltaV_arrival = 0.0;

  schedule.blocks.push_back(initblock);

  deltaV_of_schedule_init = 0.0;

  for(int i = 1; i < satnames.size(); i++){

    task_block block;

    block.arrival_time = arrival_times[i];

    block.service_duration = service_time;

    block.satname = satnames[i];

    block.departure_time = departure_times[i];

    if (deadlines.size() > 0) {

      block.demand_deadline = deadlines[i];

    }

    else {

      block.demand_deadline =  departure_times[i];
    }

    double deltaV_transfer = 0.0;

    find_optimal_trajectory_no_iter(satnames[i - 1], block.satname,
                                    departure_times[i - 1], block.arrival_time,
                                    simfile, deltaV_transfer);

    block.deltaV_arrival = deltaV_transfer;

    deltaV_of_schedule_init += deltaV_transfer;

    schedule.blocks.push_back(block);
  }

  return schedule;
}

// use local search to find the optimal_schedule

schedule_struct local_search_opt_schedule_lambert_only(double &init_deltaV, schedule_struct init_schedule, std::vector<double> dt_move,
                                                        DataFrame simfile, double service_time,std::vector<std::string> move_methods){

  double deltaVminima_so_far = init_deltaV;

  schedule_struct optimal_schedule;

  
  std::vector<double> arrival_constraints; 
  std::vector<double> departure_constraints; 


  for(int b=0 ; b < init_schedule.blocks.size(); b++){

    arrival_constraints.push_back(init_schedule.blocks[b].arrival_time);
    departure_constraints.push_back(init_schedule.blocks[b].departure_time);


  }

  while (true) {

    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;

    arma::vec deltaVs_of_neighbourhood(( init_schedule.blocks.size() * move_methods.size() * dt_move.size() ));

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
              

              if (condition_add || condition_add){
              
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

    // if no improvement is found, then stop
    if (neighbourhood_minima >= deltaVminima_so_far) {
      
      init_deltaV = deltaVminima_so_far;
      
      break;
    }

    
    //if (neighbourhood_minima == deltaVminima_so_far){
        
        
        //arma::uvec idxs = arma::regspace<arma::uvec>(0, deltaVs_of_neighbourhood.n_elem - 1);
        
        //int idxmin = static_cast<int>(arma::index_min(deltaVs_of_neighbourhood));
        //arma::uvec mask = arma::find(idxs!=idxmin); 

        //idxs = idxs.elem(mask);

        //arma::uword idx = arma::randi<arma::uword>(arma::distr_param(0, idxs.n_elem - 1));

        //int random_id = static_cast<int>(idx);

        //int worstsol = static_cast<int>(arma::index_max(deltaVs_of_neighbourhood));  


        //init_schedule = list_of_schedules[worstsol]; 
        //deltaVminima_so_far = deltaVs_of_neighbourhood[worstsol]; 

    //}  
  


    // else if improvement is found, adopt new solution
    else {

      arma::uword index_minima_uword =
          arma::index_min(deltaVs_of_neighbourhood);

      int index_minima = static_cast<int>(index_minima_uword);

      deltaVminima_so_far = std::abs(neighbourhood_minima);
      init_schedule = list_of_schedules[index_minima];

      view_schedule(list_of_schedules[index_minima]);
    }


  }//closes while

  
  optimal_schedule = init_schedule;
  return optimal_schedule;


}


schedule_struct local_search_opt_schedule(double init_deltaV,
                                          schedule_struct init_schedule,
                                          double dt_move, DataFrame simfile){

  double deltaVminima_so_far = init_deltaV;

  schedule_struct optimal_schedule;

  while (true) {

    // neighbourhood solutions
    std::vector<schedule_struct> list_of_schedules;
    arma::vec deltaVs_of_neighbourhood(init_schedule.blocks.size());
    double neighbourhood_minima;


    for (int b = 0; b < init_schedule.blocks.size(); b++) {

      schedule_struct schedule_sol = init_schedule;

      if (b == 0) {

        // move_dt(schedule_sol.blocks[b],dt_move);

        schedule_sol.blocks[b].departure_time += dt_move;
        schedule_sol.blocks[b].deltaV_arrival = 0.0;

        list_of_schedules.push_back(schedule_sol);

      }

      else {

        move_dt(schedule_sol.blocks[b], dt_move);

        schedule_sol.blocks[b].deltaV_arrival = compact_optimal_calc(
            schedule_sol.blocks[b - 1].satname, schedule_sol.blocks[b].satname,
            schedule_sol.blocks[b - 1].departure_time,
            schedule_sol.blocks[b].arrival_time, simfile);

        list_of_schedules.push_back(schedule_sol);
      }

      double totalDeltaV_of_sol = 0.0;

      for (int elem = 0; elem < schedule_sol.blocks.size(); elem++) {
        
          totalDeltaV_of_sol +=std::abs(schedule_sol.blocks[elem].deltaV_arrival);

        
      }

      deltaVs_of_neighbourhood(b) = (totalDeltaV_of_sol);

      // view_schedule(schedule_sol);
    }

    neighbourhood_minima = deltaVs_of_neighbourhood.min();

    // if no improvement is found, then stop
    if (neighbourhood_minima >= deltaVminima_so_far) {
      break;
    }

    // else if improvement is found, adopt new solution
    else {

      // std::cout<<"improvement found"<<std::endl;

      arma::uword index_minima_uword =
          arma::index_min(deltaVs_of_neighbourhood);

      int index_minima = static_cast<int>(index_minima_uword);

      deltaVminima_so_far = std::abs(neighbourhood_minima);
 

      init_schedule = list_of_schedules[index_minima];
      // view_schedule(list_of_schedules[index_minima]);
    }
  }

  optimal_schedule = init_schedule;

  return optimal_schedule;
}



void run_local_search( DataFrame simfile, std::vector<double> move_size,  
                      std::vector<std::string> moves_to_consider,
                      std::vector<std::string> sat_names_in_schedule ,
                      std::vector<double> t_depart, std::vector<double> t_arrive, 
                      double &deltaV_of_schedule, double service_time){






    //construct initial schedule using provided departure and arrival times 
    schedule_struct init_schedule = create_schedule_lambert_only(deltaV_of_schedule, t_arrive, t_depart, 
                                                                  sat_names_in_schedule,simfile, service_time);

    std::cout<< "\n"; 

    std::cout << "Initial Schedule"<< std::endl; 

    view_schedule(init_schedule);
    
    std::cout<< "\n"; 

    std::cout << "Total Delta V: "<< deltaV_of_schedule;

    // Now the optimal schedule is calculated using local search
    schedule_struct findopt_schedule = local_search_opt_schedule_lambert_only(deltaV_of_schedule, init_schedule, move_size, 
                                                                              simfile,service_time, moves_to_consider);
    
    std::cout<<"\n";

    std::cout << "Result of Local Search"<<std::endl; 

    view_schedule(findopt_schedule);

    std::cout<<"\n";

    std::cout<<"Total Delta V: "<< deltaV_of_schedule<<"\n";


    


}

































