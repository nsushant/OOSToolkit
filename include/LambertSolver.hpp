

// import libs

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>


//r1 = [r11,r12,r13], r2 = [r21,r22,r23], t and μ

arma::vec3 position_1; 
arma::vec3 position_2;
double time_of_flight;
double mu;




std::vector<arma::vec> findxy(const double &lambda, const double &T, const int &max_revolutions){

    std::vector<arma::vec> xy_pairs;

    arma::vec x_list;
    arma::vec y_list;

    ///// find xy algorithm 2 in D. Izzo "Revisiting Lambert's Problem" /////

    if (lambda > 1 ){

        std::cerr << "Error: lambda must be less than 1" << std::endl;
        return std::vector<arma::vec>();
    }

    if (T > 0.0) {

        std::cerr << "Error: T must be negative" << std::endl;
        return std::vector<arma::vec>();
    }


    double M_max = floor(T/M_PI); 
    double T_00 = acos(lambda) + lambda * sqrt(1 - lambda*lambda) + M_max * M_PI; 

    if ((T < T_00 + M_max * M_PI) && (max_revolutions > 0)){


        ///// START ALLEY ITERATIONS HERE TO FIND T_MIN /////

        double T_min; 


        if (T_min > T){

            M_max = M_max - 1;
        }


    }

    double T_1 = 2/3 * (1 - pow(lambda,3)); 

    // to be implemented 

    return xy_pairs;


}








// the lambert solver is going to return velocitiies
std::vector<arma::vec> lambert_solver(const arma::vec3& r1, const arma::vec3& r2,const double &time_of_flight, const double mu ,const int &direction, const int &max_revolutions){


    if (time_of_flight < 0){

        std::cerr << "Error: time of flight must be non-negative" << std::endl;
        return std::vector<arma::vec>();
    }

    if (mu < 0){

        std::cerr << "Error: mu must be non-negative" << std::endl;
        return std::vector<arma::vec>();
    }

    
    arma::vec3 posvec1 = r1;
    arma::vec3 posvec2 = r2;
    double tof = time_of_flight;

    // algorithm 1 D. Izzo "Revisiting Lambert's Problem" 
    arma::vec3 c = posvec2 - posvec1;

    // vec norms
    double c_norm = arma::norm(c);
    double r1_norm = arma::norm(posvec1);
    double r2_norm = arma::norm(posvec2);


    //calc s 
    double s = 0.5 * (r1_norm + r2_norm + c_norm);
    

    arma::vec3 dirvec_r1= posvec1 / r1_norm; // unit vector in the direction of r1
    arma::vec3 dirvec_r2 = posvec2 / r2_norm; // unit vector in the direction of r2
    arma::vec h_vec_dir = arma::cross(dirvec_r1, dirvec_r2); // direction of the angular momentum vector

    double lambda = sqrt(1 - c_norm / s);


    arma::vec dirvec_t1;
    arma::vec dirvec_t2;

    // direction based change in lambda
    if (r1(0)*r2(1) - r1(1)*r2(0) < 0){
        lambda = -lambda;

        arma::vec dirvec_t1 = arma::cross(dirvec_r1,h_vec_dir); 
        arma::vec dirvec_t2 = arma::cross(dirvec_r2,dirvec_r2); 

    }

    else
    {
        arma::vec dirvec_t1 = arma::cross(h_vec_dir, dirvec_r1); 
        arma::vec dirvec_t2 = arma::cross(h_vec_dir,h_vec_dir); 

    }

    double T = sqrt( 2*mu / pow(s,3) ) * tof;


    ///// outputs of findxy ////// <- needs to be implemented
    //findxy();
    
    // assume equal length xlist and ylist
    arma::vec xlist; 
    arma::vec ylist; 


    double gamma = sqrt(mu * s / 2); 
    double rho = (r1_norm - r2_norm)/ c_norm ; 
    double sigma = sqrt(1 - pow(rho,2)); 

    std::vector<arma::vec> velocity_solutions;

    
    // looping over all xy pairs returned from findxy
    for (int elem = 0; elem < xlist.n_elem; elem++){

        double x_loop = xlist(elem);
        double y_loop = ylist(elem);


        double v_r1 = gamma*( (lambda*y_loop - x_loop) -  rho*( lambda*y_loop + x_loop) )  / r1_norm;
        double v_r2 = -gamma*( (lambda*y_loop - x_loop) +  rho*( lambda*y_loop + x_loop) )  / r2_norm;

        double v_t1 = gamma*sigma*( y_loop + lambda*x_loop ) / r1_norm;
        double v_t2 = -gamma*sigma*(y_loop + lambda*x_loop) / r2_norm;

        arma::vec v1 = v_r1*dirvec_r1 + v_t1*dirvec_t1;
        arma::vec v2 = v_r2*dirvec_r2 + v_t2*dirvec_t2;

        velocity_solutions.push_back(arma::join_cols(v1,v2));

    }

    return velocity_solutions;


}