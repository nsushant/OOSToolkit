#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>

#include "Local_search.hpp"
#include "data_access_lib.hpp"
#include "LambertSolver.hpp"
#include "Nbody.hpp"

/*
 The outputs of local search algorithms need to be compared with those of exact methods, so as to 
 quantify the difference between the solution obtained and the optimal solution (to calculate the optimality gap).  

author: Sushanta Nigudkar
date: 11/2025

*/


// The aim is to use dynamic programming to find an optimal schedule. 


bool is_feasible_sol(task_block block1, task_block block2, DataFrame simfile, double t1, double t2){

    

    std::string name1 = block1.satname; 
    std::string name2 = block2.satname; 


    std::vector<std::string> satnames = simfile["name"];

    arma::uvec sat_idxs_1= find_idxs_of_match(satnames,name1);
    arma::uvec sat_idxs_2 = find_idxs_of_match(satnames,name2);


    arma::vec t_sat = simfile.getNumeric("time_s");
    
    arma::vec ts_sat1 = t_sat.elem(sat_idxs_1);
    arma::vec ts_sat2 = t_sat.elem(sat_idxs_2);
    

    arma::uvec sat1_t1 = arma::find(ts_sat1 == t1);
    arma::uvec sat2_t2 = arma::find(ts_sat2 == t2);  

    arma::vec x_sat = simfile.getNumeric("x");
    arma::vec y_sat = simfile.getNumeric("y");
    arma::vec z_sat = simfile.getNumeric("z");
    

    arma::vec tmp1 = x_sat.elem(sat1_t1);
    arma::vec tmp2 = y_sat.elem(sat1_t1);
    arma::vec tmp3 = z_sat.elem(sat1_t1);

    arma::vec3 vec_r1 = {tmp1(0), tmp2(0), tmp3(0)};

    arma::vec tmp4 = x_sat.elem(sat2_t2);
    arma::vec tmp5 = y_sat.elem(sat2_t2);
    arma::vec tmp6 = z_sat.elem(sat2_t2);

    arma::vec3 vec_r2 = {tmp4(0),tmp5(0),tmp6(0)};
    
    double r1_norm = arma::norm(vec_r1); 
    double r2_norm = arma::norm(vec_r2);


    // algorithm 1 D. Izzo "Revisiting Lambert's Problem" 
    arma::vec3 vec_c = vec_r2 - vec_r1;

    // vec norms
    double c = arma::norm(vec_c);
    
    double s = (r1_norm + r2_norm + c) / 2.0;

    double lambda = sqrt(1.0 - c / s);

    double tof = t2 - t1; 

    double T = sqrt( 2.0*MU_EARTH / (s*s*s) ) * tof;
    
    double T_1 = 2.0/3.0 * (1.0 - pow(lambda,3));
    
    if(T <= T_1){

        return false;
    }
    
    else{

        return true;

    }
    
}



void finding_individual_minimas_dynamic_programming(schedule_struct &init_schedule, DataFrame Simfile, double dt){

    
    // tabulate this table with solutions as they are calculated 
    
    // table[arrival_time][departure_time]
    // lets discretize into dt second intervals 
    // 
    // if we start with a schedule that lasts 60000 seconds we will get 6000 timesteps
    // out of these 6000 if dt = 100. 
    // 
    // There will be further reductions due to the time occupied by servicing tasks themselves 
    // (times at which we can neither depart nor arrive)
    // 
    // Additionally we can check if a solution has a time of flight = T1 and omit these solutions 

    // The table will have a maximum of 4 cols (arrival time, departure time, name, delta V)
    // and a maximum of arrival_time_of_last_block/100 rows.


    int numtimes = std::ceil(init_schedule.blocks[(int)init_schedule.blocks.size() - 1].arrival_time/dt); 
    int numrows = numtimes * numtimes; 

    //view_schedule(init_schedule); 

    //std::cout << "dt : " << dt << std::endl;
    std::cout << "numrows : " << numrows << std::endl;

    arma::mat table_of_sols(numrows,4,arma::fill::zeros);

    schedule_struct optimal_schedule;
    
    double sols = 0; 

    for(int b = 1 ; b < init_schedule.blocks.size(); b++){

        task_block block1 = init_schedule.blocks[b-1]; 
        task_block block2 = init_schedule.blocks[b];

        double d_constraint = block1.departure_time; 
        double a_constraint = block2.arrival_time; 
        
            
        task_block optblock;
        optblock = block2; 
        

        task_block initblock;
        initblock = block1; 
        
        std::vector<double> dopt; 
        std::vector<double> aopt; 
        std::vector<double> deltaopt; 
        

        for(double a = a_constraint; a > d_constraint; a-=dt){
            
            for(double d = d_constraint; d <= a; d+=dt){
            
            //std::cout << "d,a : "<< d <<" , "<<a << std::endl;  

            if(is_feasible_sol(block1, block2, Simfile, d, a)){ 


                arma::rowvec current_prob = {a,d,b}; 

                double tol = 1e-9; 


                //arma::uvec idxs_sat = arma::find(arma::all( arma::abs(table_of_sols.cols(0,2).each_row() - current_prob) < tol, 1 )); 

                arma::mat matches_with_current_prob = table_of_sols.cols(0,2).each_row() - current_prob; 
                arma::vec sum_mat_cols = arma::sum(matches_with_current_prob,1);
                arma::uvec idxs_sat = arma::find(sum_mat_cols == 0); 
                
                

                double DeltaVsol = 0; 

                if ((int)idxs_sat.n_elem == 0){
            
                    //std::cout << "computing new soln" << "\n"; 

                    find_optimal_trajectory_no_iter(block1.satname, block2.satname, d, a, Simfile, DeltaVsol); 
            
                    //check if solution exists in lookup table

                    table_of_sols(sols,0) = a; 
                    table_of_sols(sols,1) = d; 
                    table_of_sols(sols,3) = DeltaVsol; 
                    table_of_sols(sols,2) = b; 

                    //std::cout <<"wrote solution"<<std::endl;
                    sols+=1;
                }

                else{
                    
                    //std::cout << "trying to read precomputed sol" << "\n";

                    arma::mat readsol = table_of_sols.rows(idxs_sat);

                    DeltaVsol = readsol(0,3); 

                }

                
                deltaopt.push_back(DeltaVsol); 
                dopt.push_back(d);
                aopt.push_back(a); 

                }
            
            }

        }



        // find time pair that gives minimum delta V 
        std::vector<double>::iterator it = std::min_element(deltaopt.begin(), deltaopt.end());
        // Get index
        std::size_t idxmin = std::distance(deltaopt.begin(), it);

        optblock.deltaV_arrival = deltaopt[idxmin];
        optblock.arrival_time = aopt[idxmin]; 
        initblock.departure_time = dopt[idxmin]; 


        if (b==1){

            optimal_schedule.blocks.push_back(initblock); 
        
        }


        else{

            optimal_schedule.blocks[ optimal_schedule.blocks.size() - 1 ].departure_time = dopt[idxmin]; 
            
        }
        
        optimal_schedule.blocks.push_back(optblock);

    }
                
        
    init_schedule = optimal_schedule;


} 

