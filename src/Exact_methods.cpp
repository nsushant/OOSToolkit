#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

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



/*
double calc_task_cost(int n,DataFrame& Simfile,double&DeltaVsol){

    
    task_block block1 = init_schedule.blocks[n-1]; 
    task_block block2 = init_schedule.blocks[n];

    find_optimal_trajectory_no_iter(block1.satname, block2.satname, block1.departure_time, block2.arrival_time, Simfile, DeltaVsol); 
    

    calc_task_cost(n-1, Simfile, DeltaVsol); 

}
*/









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


    //int numtimes = std::ceil(init_schedule.blocks[(int)init_schedule.blocks.size() - 1].arrival_time/dt); 
    //int numrows = numtimes * numtimes; 

    //view_schedule(init_schedule); 

    //std::cout << "dt : " << dt << std::endl;
    //std::cout << "numrows : " << numrows << std::endl;

    //arma::mat table_of_sols(numrows,4,arma::fill::zeros);
    
    // first lets assume the following constraints apply 
    // d_i < a_{i+1}
    // a_i < d_i 
    // d_i > a_i + s_i 
    //
    //
    // If we now consider fixed arrival times b_i to get additional constraints 
    // 
    // a_i < b_i 
    // d_i > b_i + s_i  
    //
    // for each solution we would need the block number, deprature time, arrival time, delta V 
    //  
    // suppose we get 
    //
    //
    //Table of sols [a, d, b] = delta V
    //
    std::map<std::vector<double>,double> table_of_sols; 
	
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
            
            for(double d = d_constraint; d < a; d+=dt){
            
            //std::cout << "d,a : "<< d <<" , "<<a << std::endl;  

            if(is_feasible_sol(block1, block2, Simfile, d, a)){ 


		        std::vector<double> current_prob = {a,d,(double)b}; 
                
                                   
                bool match_found = (table_of_sols.find(current_prob) != table_of_sols.end());  
                
                double DeltaVsol = 0; 

                if (match_found == false){
            
                    //std::cout << "computing new soln" << "\n"; 

                    find_optimal_trajectory_no_iter(block1.satname, block2.satname, d, a, Simfile, DeltaVsol); 
            
                    //check if solution exists in lookup table
                        
                    table_of_sols[current_prob] = DeltaVsol; 

                    //std::cout <<"wrote solution"<<std::endl;
                    sols+=1;
                }

                else{
                    
                    //std::cout << "trying to read precomputed sol" << "\n";

                    DeltaVsol = table_of_sols[current_prob]; 

                }

                
                deltaopt.push_back(DeltaVsol); 
                dopt.push_back(d);
                aopt.push_back(a); 

            }
            
          }

        }


        // find time pair that gives minimum delta V 
        auto it = std::min_element(deltaopt.begin(), deltaopt.end());
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





double calculate_deltav_upto_thisblock(int n, schedule_struct&sched_in, DataFrame simfile, std::map<std::vector<double>, double>&lookup){


    if(n==0){

        return 0.0; 
    }


    std::vector<double> key = {(double)n,sched_in.blocks[n-1].departure_time,sched_in.blocks[n].arrival_time};

    
    double deltaV_current;

    if (lookup.find(key) != lookup.end()){

                
        deltaV_current+=lookup[key];

        sched_in.blocks[n].deltaV_arrival = deltaV_current; 

    }

    else{
            
        // if you want to include considered wait-times uncomment the line below

        //deltaV_current += compact_optimal_calc( sched_in.blocks[n-1],sched_in.blocks[n],
            //                                      sched_in[n-1].departure_time, sched_in[n].arrival_time,simfile);
         
        
        //std::cout <<"n,d,a"<< n <<","<<sched_in.blocks[n-1].departure_time<<","<< sched_in.blocks[n].arrival_time << "\n";  

        find_optimal_trajectory_no_iter(sched_in.blocks[n-1].satname, sched_in.blocks[n].satname, 
                                        sched_in.blocks[n-1].departure_time, sched_in.blocks[n].arrival_time, 
                                        simfile, deltaV_current);

        sched_in.blocks[n].deltaV_arrival = deltaV_current; 


            
        lookup[key] = deltaV_current; 

    }    


    return deltaV_current += calculate_deltav_upto_thisblock(n-1, sched_in,simfile,lookup);       
    
    
}



int dynamic_program_fixed_tasksize_Tfixed(schedule_struct& sched_init, double dt,int b_reached,int b_lim,double Tmax,DataFrame simfile){




    if (b_reached == b_lim){
            
        std::cout << " dynamic program terminated" << std::endl;  
        return 0.0; 
    
    }

    // 
    // constraints 
    // a_i < d_{i-1} 
    // a_i = d_i - s_i 
    // a_i > 0 
    // s_i > 0
    // d_i     

    std::map<std::vector<double>, double>lookup;
   

        
    double deltav_so_far = 0.0; 

    for(int itb = 1 ; itb <= b_lim; itb++)
    {

        deltav_so_far += sched_init.blocks[itb].deltaV_arrival; 

    }


    for(int b = b_reached ; b>=1 ; b--)
    {

        std::cout << "b_reached : " << b << std::endl; 
        double dep_init = (b > 1) ? sched_init.blocks[b-2].departure_time : 0.0; 
    
        double arr_init = sched_init.blocks[b].arrival_time;  

        schedule_struct sched_sub = sched_init; 
        
        if (b > 1){

        while( sched_sub.blocks[b-1].arrival_time > dep_init )
        {
            
            sched_sub.blocks[b-1].departure_time-=dt ;
            
            if( (b-1) > 0 )
            {
                sched_sub.blocks[b-1].arrival_time-=dt;   
               
                if(sched_sub.blocks[b-2].departure_time >= sched_sub.blocks[b-1].arrival_time)
                {
                    break;

                }

            } 
            

            bool feasibility = is_feasible_sol(sched_sub.blocks[b-1], sched_sub.blocks[b], simfile, 
                                sched_sub.blocks[b-1].departure_time, sched_sub.blocks[b].arrival_time);
            
            if(feasibility == false){

                continue; 
            }


            double deltav_upto_block = calculate_deltav_upto_thisblock(b_reached, sched_sub, simfile, lookup);
                
                       
            if (deltav_upto_block < deltav_so_far)
            {           
                sched_init = sched_sub;
                deltav_so_far = deltav_upto_block;   

            }
                
                

        }

        }
        
        
        

        schedule_struct sched_add = sched_init; 
        std::cout <<"arr_int"<<arr_init<<std::endl;
        while( sched_add.blocks[b-1].departure_time < arr_init)
        {

            sched_add.blocks[b-1].departure_time +=dt;
            sched_add.blocks[b-1].arrival_time   +=dt;
        

            if(sched_add.blocks[b].arrival_time <= sched_add.blocks[b-1].departure_time){
                break; 
            }
             


            
            bool feasibility = is_feasible_sol(sched_add.blocks[b-1], sched_add.blocks[b], simfile, 
                                sched_add.blocks[b-1].departure_time, sched_add.blocks[b].arrival_time);
            
            if(feasibility == false){

                continue; 
            }

            
            double deltav_upto_block = calculate_deltav_upto_thisblock(b_reached, sched_add, simfile, lookup);

            
            if (deltav_upto_block < deltav_so_far){

                //sched_init.blocks[b-1].departure_time = sched_add.blocks[b-1].departure_time; 
                //sched_init.blocks[b-1].arrival_time = sched_add.blocks[b-1].arrival_time; 
                sched_init = sched_add; 
                deltav_so_far = deltav_upto_block; 
            }
            
        

        }

    }
     
    
    b_reached +=1; 

    return dynamic_program_fixed_tasksize_Tfixed(sched_init,dt,b_reached,b_lim,Tmax,simfile);
  

}

















