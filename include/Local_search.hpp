#pragma once 
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include "Trajectory_selection.hpp"
#include "data_access_lib.hpp"

/*
Author : S. Nigudkar (2025)

Functions to initialize schedules and find the optimal schedule using local search. 

*/

struct task_block{

    double arrival_time; 
    double departure_time;
    std::string satname; 
    double deltaV_arrival; 


};

struct schedule_struct{


    std::vector<task_block> blocks;


};

void move_dt(task_block &tb, double dt){

    tb.arrival_time -= dt; 
    tb.departure_time -= dt;

}

void move_dt2(std::vector<task_block> &blocks, int b_index, double dt){

    if (b_index == 0){
        
        blocks[b_index].departure_time += dt/2 - ((int)dt/2 %100);

    }


    else{

        if (blocks[b_index].arrival_time > blocks[b_index-1].departure_time){

                blocks[b_index].arrival_time -= dt; 
                //blocks[b_index].arrival_time += dt; 

            
            if (blocks[b_index].departure_time == 0){
                
                blocks[b_index].departure_time += 0;

            }
            
            else{
                
                blocks[b_index].departure_time += dt/2 - ((int)dt/2 %100) ;
           
            }

        
        }


    }


}



void move_add_arrival(std::vector<task_block> &blocks, int b_index, double dt){

    if (b_index < (blocks.size() - 1)){

    if (blocks[b_index].arrival_time + dt < blocks[b_index].departure_time){
        
        blocks[b_index].arrival_time += dt;

    } 
    
    }

    else{

        blocks[b_index].arrival_time += dt;

    }    

}
void move_add_departure(std::vector<task_block> &blocks, int b_index, double dt){

    if (b_index + 1 <= blocks.size()-1){

    if (blocks[b_index].departure_time + dt < blocks[b_index + 1].arrival_time){
        
        blocks[b_index].departure_time += dt;  

    }
    
    }

    else{

        blocks[b_index].departure_time += 0; 

    }

}

void move_sub_departure(std::vector<task_block> &blocks, int b_index, double dt){ 

    if (b_index < (blocks.size()-1)){

        if (blocks[b_index].departure_time - dt > blocks[b_index].arrival_time){

            blocks[b_index].departure_time -= dt;
        }




    }

}
void move_sub_arrival(std::vector<task_block> &blocks, int b_index, double dt){


    if ((b_index - 1) < (blocks.size()-1)){

        if (blocks[b_index].arrival_time - dt > blocks[b_index - 1].departure_time){

            blocks[b_index].arrival_time -= dt;
        }

    }


    else
    {   
        if (blocks[b_index].arrival_time - dt > 0){

        blocks[b_index].arrival_time -= dt;

        }

    }
    
}


void move_wrapper(std::vector<task_block> &blocks, int b_index, double dt, std::string method){

    if (method == "add arrival"){ move_add_arrival(blocks, b_index, dt);}
    if (method == "add departure"){ move_add_departure(blocks, b_index, dt);}
    if (method == "sub arrival"){ move_sub_arrival(blocks, b_index, dt);}
    if (method == "sub departure"){ move_sub_departure(blocks, b_index, dt);}
    if (method == "move_dt2"){move_dt2(blocks, b_index, dt);}

}



void view_schedule(schedule_struct schedule_to_print){

    // Define column widths
    const int w1 = 15;   // Task
    const int w2 = 15;  // Arrival
    const int w3 = 15;  // Departure
    const int w4 = 13;  // Fuel cost

    // Print header line
    std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";
    std::cout << "| " << std::left << std::setw(w1 - 1) << "Task"
              << "| " << std::setw(w2 - 1) << "Arrival Time"
              << "| " << std::setw(w3 - 1) << "Departure Time"
              << "| " << std::setw(w4 - 1) << "Delta V"
              << "|\n";
    std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";

    // Print each task
    for (const auto& t : schedule_to_print.blocks) {
        std::cout << "| " << std::left << std::setw(w1 - 1) << t.satname
                  << "| " << std::setw(w2 - 1) << t.arrival_time
                  << "| " << std::setw(w3 - 1) << t.departure_time
                  << "| " << std::right << std::setw(w4 - 3)
                  << std::fixed << std::setprecision(2) << t.deltaV_arrival
                  << "  |\n";
    }

    std::cout << std::string(w1 + w2 + w3 + w4 + 7, '-') << "\n";


}

schedule_struct create_schedule(double &deltaV_of_schedule_init, std::vector<double> arrival_times,std::vector<double> departure_times, std::vector<std::string> satnames,DataFrame simfile){

    // initial block for the service depot
    
    schedule_struct schedule; 

    task_block initblock; 

    initblock.arrival_time = arrival_times[0]; 

    initblock.departure_time = departure_times[0]; 

    initblock.satname = satnames[0];

    initblock.deltaV_arrival = 0.0; 

    schedule.blocks.push_back(initblock); 

    deltaV_of_schedule_init = 0.0; 

    for (int i = 1 ; i < satnames.size() ; i++){

        task_block block; 
        
        block.arrival_time = arrival_times[i];

        block.departure_time = departure_times[i]; 

        block.satname = satnames[i];

        double deltaV_transfer = compact_optimal_calc(satnames[i-1], block.satname, departure_times[i-1], block.arrival_time, simfile);

        block.deltaV_arrival = deltaV_transfer; 

        deltaV_of_schedule_init+= deltaV_transfer; 

        schedule.blocks.push_back(block);

    }


    return schedule;


}






schedule_struct create_schedule_lambert_only(double &deltaV_of_schedule_init, std::vector<double> arrival_times,std::vector<double> departure_times, std::vector<std::string> satnames,DataFrame simfile, double service_time){

    // initial block for the service depot
    
    schedule_struct schedule; 

    task_block initblock; 

    initblock.arrival_time = arrival_times[0]; 

    initblock.departure_time = departure_times[0]; 

    initblock.satname = satnames[0];

    initblock.deltaV_arrival = 0.0; 

    schedule.blocks.push_back(initblock); 

    deltaV_of_schedule_init = 0.0; 

    for (int i = 1 ; i < satnames.size() ; i++){

        task_block block; 
        
        block.arrival_time = arrival_times[i];

        block.departure_time = departure_times[i]; 

        //block.departure_time = block.arrival_time + service_time; 

        block.satname = satnames[i];

        double deltaV_transfer = 0.0;

        find_optimal_trajectory_no_iter(satnames[i-1], block.satname, departure_times[i-1], block.arrival_time, simfile, deltaV_transfer);         

        block.deltaV_arrival = deltaV_transfer; 

        deltaV_of_schedule_init+= deltaV_transfer; 

        schedule.blocks.push_back(block);

    }

    return schedule;

}

schedule_struct local_search_opt_schedule_lambert_only(double &init_deltaV, schedule_struct init_schedule, double dt_move, DataFrame simfile, double service_time ,std::vector<std::string> move_methods){

    double deltaVminima_so_far = init_deltaV; 

    schedule_struct optimal_schedule; 

    while (true){
    
        //neighbourhood solutions
        std::vector<schedule_struct> list_of_schedules;
        arma::vec deltaVs_of_neighbourhood((init_schedule.blocks.size()*move_methods.size())); 

        double neighbourhood_minima; 

        int solnum_neighbourhood = 0; 

        for (int m = 0 ; m < move_methods.size() ; m++ ){

        for (int b = 0 ; b < init_schedule.blocks.size() ; b++ ){


                schedule_struct schedule_sol = init_schedule;

                move_wrapper(schedule_sol.blocks, b, dt_move,move_methods[m]);

                if (b == 0){

                    //move_dt(schedule_sol.blocks[b],dt_move); 
                    //move_dt2(schedule_sol.blocks, b, dt_move); 

                    schedule_sol.blocks[b].deltaV_arrival = 0.0; 

                    list_of_schedules.push_back(schedule_sol); 


                }

                
                else {

                //move_dt(schedule_sol.blocks[b],dt_move);
                //move_dt2(schedule_sol.blocks, b, dt_move); 

                double DeltaVMinimaopt; 


                find_optimal_trajectory_no_iter(schedule_sol.blocks[b-1].satname, schedule_sol.blocks[b].satname, schedule_sol.blocks[b-1].departure_time, schedule_sol.blocks[b].arrival_time, simfile, DeltaVMinimaopt); 
                
                //find_optimal_trajectory_no_iter(schedule_sol.blocks[b-1].satname, schedule_sol.blocks[b].satname, schedule_sol.blocks[b-1].arrival_time + service_time, schedule_sol.blocks[b].arrival_time, simfile, DeltaVMinimaopt); 

                schedule_sol.blocks[b].deltaV_arrival = DeltaVMinimaopt;
                
                list_of_schedules.push_back(schedule_sol); 
            
                }


                double totalDeltaV_of_sol = 0.0; 

                for (int elem = 0 ; elem < schedule_sol.blocks.size() ; elem++){

                    
                    totalDeltaV_of_sol += std::abs(schedule_sol.blocks[elem].deltaV_arrival); 

                }


                deltaVs_of_neighbourhood(solnum_neighbourhood) = (totalDeltaV_of_sol);

                if (totalDeltaV_of_sol == 0){ std::cout<<"0 delta V";}

                solnum_neighbourhood+=1;

                //view_schedule(schedule_sol);

        }

        }



        neighbourhood_minima = deltaVs_of_neighbourhood.min(); 
        //std::cout<< "neighbourhood_minima : "<< neighbourhood_minima << std::endl;
        //std::cout<< "neighbourhood_maxima : "<< deltaVs_of_neighbourhood.max() << std::endl;
        //std::cout<< "size : " << deltaVs_of_neighbourhood.size()<< std::endl;

        // if no improvement is found, then stop
        if (neighbourhood_minima >= deltaVminima_so_far){
            init_deltaV = deltaVminima_so_far; 
            break;
        }
        
        //else if improvement is found, adopt new solution 
        else{

            arma::uword index_minima_uword = arma::index_min(deltaVs_of_neighbourhood); 

            int index_minima = static_cast<int>(index_minima_uword);        

            deltaVminima_so_far = std::abs(neighbourhood_minima);
            init_schedule = list_of_schedules[index_minima];    
            
            //view_schedule(list_of_schedules[index_minima]); 


        }

    }

    optimal_schedule = init_schedule; 

    return optimal_schedule;

}









schedule_struct local_search_opt_schedule(double init_deltaV, schedule_struct init_schedule, double dt_move, DataFrame simfile){


    double deltaVminima_so_far = init_deltaV; 

    schedule_struct optimal_schedule; 

    while (true){
    
        //neighbourhood solutions
        std::vector<schedule_struct> list_of_schedules;
        arma::vec deltaVs_of_neighbourhood(init_schedule.blocks.size()); 
        double neighbourhood_minima; 

        for (int b = 0 ; b < init_schedule.blocks.size() ; b++ ){

                schedule_struct schedule_sol = init_schedule;
                
                if (b == 0){

                    //move_dt(schedule_sol.blocks[b],dt_move); 

                    schedule_sol.blocks[b].departure_time += dt_move;
                    schedule_sol.blocks[b].deltaV_arrival = 0.0; 

                    list_of_schedules.push_back(schedule_sol); 

                }

                else {

                move_dt(schedule_sol.blocks[b],dt_move);

                schedule_sol.blocks[b].deltaV_arrival =  compact_optimal_calc(schedule_sol.blocks[b-1].satname, schedule_sol.blocks[b].satname, schedule_sol.blocks[b-1].departure_time, schedule_sol.blocks[b].arrival_time, simfile);

                list_of_schedules.push_back(schedule_sol); 
            
        
                }

                double totalDeltaV_of_sol = 0.0; 


                for (int elem = 0 ; elem < schedule_sol.blocks.size() ; elem++){

                    
                    totalDeltaV_of_sol += std::abs(schedule_sol.blocks[elem].deltaV_arrival); 

                }

                deltaVs_of_neighbourhood(b) = (totalDeltaV_of_sol);

                //view_schedule(schedule_sol);

        }


        neighbourhood_minima = deltaVs_of_neighbourhood.min(); 

        // if no improvement is found, then stop
        if (neighbourhood_minima >= deltaVminima_so_far){
            break;
        }
        
        //else if improvement is found, adopt new solution 
        else{

            //std::cout<<"improvement found"<<std::endl;
            

            arma::uword index_minima_uword = arma::index_min(deltaVs_of_neighbourhood); 
        
            int index_minima = static_cast<int>(index_minima_uword);        

            deltaVminima_so_far = std::abs(neighbourhood_minima);
            init_schedule = list_of_schedules[index_minima];    

            //view_schedule(list_of_schedules[index_minima]); 


        }

    }

    optimal_schedule = init_schedule; 

    return optimal_schedule;

}






