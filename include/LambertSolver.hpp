

/*

The lambert solver has been adopted from pykep with parts directly obtained from both the 
github repo located at (https://github.com/esa/pykep/blob/832db1aecd2cf1b6a752f86e50daef3c02ffbdcb/src/lambert_problem.cpp). 

The original publication from D.IZZO detailing the algorithm used to solve lambert's probelm is available 
at (https://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2014-RevisitingLambertProblem.pdf) with a deviation 
in algorithm 1 noted between the author's github repo and the published algorithm. the correction has been applied in
the below implementation ( for details on this see lines 241-264 below and the equations of i_t2 in the springer 
publication's listing of Algorithm 1.) 


Deviations from pykep 

1. The use of armadillo for linear algebra and vector calculations.  
2. Replacement of boost functions with appropriate equivalents. 
3. Removal of dependencies on pre-written vector calculation headers in pykep
4. Addition of minimal error handling lines. 

*/


#pragma once

// import libs

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream> 
#include <armadillo>
#include <algorithm>


//r1 = [r11,r12,r13], r2 = [r21,r22,r23], t and μ

arma::vec3 position_1; 
arma::vec3 position_2;
double time_of_flight;
double mu;


void taking_derivatives( double &DT,  double &DDT, double &DDDT, double const lambda, const double x, const double T){


    // safe guards to avoid sqrt of negative or division by zero
    //const double tiny = 1e-14;
    double l2 = lambda * lambda;
    double l3 = l2 * lambda;

    // denominator 1 - x^2 must not be zero
    double umx2 = 1.0 - x * x;
    double y = sqrt(1.0 - l2 * umx2);

    double y2 = y * y;
    double y3 = y2 * y;
    
    DT = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y);
    DDT = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3);
    DDDT = 1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2);


}



void tof_lagrange(double &tof, const double x, const int N, const double lambda)
{
    // directly from pykep

    // here 'x' is the universal veriable 
    double a = 1.0 / (1.0 - x * x);

    if (a > 0) 
    {
        double alfa = 2.0 * acos(x);
        double beta = 2.0 * asin(sqrt(lambda*lambda / a));
        if (lambda < 0.0) beta = -beta;

        tof = ((a * sqrt(a) * ((alfa - sin(alfa)) - (beta - sin(beta)) + 2.0 * M_PI * N)) / 2.0);

    } else {
        double alfa = 2.0 * acosh(x);
        double beta = 2.0 * asinh(sqrt(-(lambda*lambda) / a));

        if (lambda < 0.0) beta = -beta;

        tof = (-a * sqrt(-a) * ((beta - sinh(beta)) - (alfa - sinh(alfa))) / 2.0);
    }
}




double hypergeometricF(double z, double tol)
{

    // directly from pykep

    double Sj = 1.0;
    double Cj = 1.0;
    double err = 1.0;
    double Cj1 = 0.0;
    double Sj1 = 0.0;
    int j = 0;

    
    while (err > tol) {
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
        Sj1 = Sj + Cj1;
        err = fabs(Cj1);
        Sj = Sj1;
        Cj = Cj1;
        j = j + 1;
    }

    return Sj;

}






void calc_time_of_flight(double &time_of_flight,const double x, const int N , const double lambda){

    // directly adapted from pykep
    double battin = 0.01;
    double lagrange = 0.2;
    double dist = fabs(x - 1);

    if (dist < lagrange && dist > battin) { // We use Lagrange tof expression
        tof_lagrange(time_of_flight, x, N,lambda);
        return;
    }
    double K = lambda*lambda;
    double E = x * x - 1.0;
    double rho = fabs(E);
    double z = sqrt(1 + K * E);
    if (dist < battin) { // We use Battin series tof expression
        double eta = z - lambda * x;
        double S1 = 0.5 * (1.0 - lambda - x * eta);
        double Q = hypergeometricF(S1, 1e-11);
        Q = 4.0 / 3.0 * Q;
        time_of_flight = (eta * eta * eta * Q + 4.0 * lambda * eta) / 2.0 + N * M_PI / pow(rho, 1.5);
        return;
    } else { // We use Lancaster tof expresion
        double y = sqrt(rho);
        double g = x * z - lambda * E;
        double d = 0.0;
        if (E < 0) {
            double l = acos(g);
            d = N * M_PI + l;
        } else {
            double f = y * (z - lambda * x);
            d = log(f + g);
        }
        time_of_flight = (x - lambda * z - d / y) / E;
        return;
    }
}


int iterate_householder(const double T, double &x0, const int N, const double eps, const int iter_max,const double lambda)
{

    int it = 0;
    double err = 1.0;
    double xnew = 0.0;
    double tof = 0.0, delta = 0.0, DT = 0.0, DDT = 0.0, DDDT = 0.0;
    double x_orig = x0;
    while ((err > eps) && (it < iter_max)) {
        calc_time_of_flight(tof, x0, N, lambda);
        taking_derivatives(DT, DDT, DDDT, lambda, x0, tof);
        delta = tof - T;
        double DT2 = DT * DT;
        xnew = x0 - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0);
        err = fabs(x0 - xnew);
        x0 = xnew;
        it++;
    }
    return it;
}



// the lambert solver is going to return velocitiies v1, v2 which have to be applied at t1,t2 respectively 
// inputs - r1 (t1) , r2 (t2), time of flight (t2-t1), mu, direction (1 for prograde, -1 for retrograde), max_revolutions (default 0)
std::vector<arma::vec> lambert_solver(const arma::vec3 r1, const arma::vec3 r2,const double time_of_flight, const double mu ,const int retrograde, const int max_revolutions){

    if (time_of_flight <=0){

        std::cerr << "Error: time of flight must be non-negative and non-zero" << std::endl;
        return std::vector<arma::vec>();
    }

    if (mu <= 0){

        std::cerr << "Error: mu must be non-negative and non-zero" << std::endl;
        return std::vector<arma::vec>();
    }


    arma::vec3 vec_r1 = r1;
    arma::vec3 vec_r2 = r2;
    double tof = time_of_flight;
    int max_revs = max_revolutions;
    double direction = retrograde;

    // algorithm 1 D. Izzo "Revisiting Lambert's Problem" 
    arma::vec3 vec_c = vec_r2 - vec_r1;

    // vec norms
    double c = arma::norm(vec_c);
    double r1_norm = arma::norm(vec_r1);
    double r2_norm = arma::norm(vec_r2);


    //calc s 
    double s = (r1_norm + r2_norm + c) / 2.0 ;
    

    arma::vec3 dirvec_r1= vec_r1 / r1_norm; // unit vector in the direction of r1
    arma::vec3 dirvec_r2 = vec_r2 / r2_norm; // unit vector in the direction of r2
    
    arma::vec h_vec_dir = arma::cross(dirvec_r1, dirvec_r2); // direction of the angular momentum vector
    h_vec_dir = h_vec_dir / arma::norm(h_vec_dir); // unit vector in the direction of h


    double lambda = sqrt(1.0 - c / s);

    if (h_vec_dir(2) == 0){

        std::cerr << "The angular momentum has no z component, cannot determine direction of orbit (cw or anti-cw)" << std::endl;

    }


    arma::vec dirvec_t1;
    arma::vec dirvec_t2;

    if (h_vec_dir(2) < 0){

        // motion is clockwise
        // long way solution 
        // transfer angle is greater than pi 

        lambda = -lambda;

        arma::vec t1_cross = arma::cross(dirvec_r1,h_vec_dir);
        dirvec_t1 = t1_cross / arma::norm(t1_cross); // unit vector in the direction of t1
       
        arma::vec t2_cross = arma::cross(dirvec_r2,h_vec_dir); 
        dirvec_t2 = t2_cross / arma::norm(t2_cross); // unit vector in the direction of t1

    }

    else
    {
        arma::vec t1_cross= arma::cross(h_vec_dir, dirvec_r1); 
        dirvec_t1 = t1_cross / arma::norm(t1_cross); // unit vector in the direction of t1

        arma::vec t2_cross = arma::cross(h_vec_dir,dirvec_r2); 
        dirvec_t2 = t2_cross / arma::norm(t2_cross); // unit vector in the direction of t1

    }

    
    if (direction == 1){

        // if cw direction is imposed

        lambda = -1.0 * lambda;
        dirvec_t1 = -1.0 * dirvec_t1;
        dirvec_t2 = -1.0 * dirvec_t2;

    }

    double lambda_squared = lambda*lambda;
    double lambda_cubed = lambda_squared * lambda;


    // calculating non dimensional time of flight T
    double T = sqrt( 2.0*mu / (s*s*s) ) * tof;

    /*

    std::cout << "Non-dimensional time of flight T : " << T << std::endl;
    std::cout << "lambda : " << lambda << std::endl;
    std::cout << "chord : " << c << std::endl;
    std::cout << "s : " << s << std::endl;

    std::cout << "finding all 'x' for given lambda, T " << std::endl;
     */

    // maximum number of possible revs given travel time
    int M_max = static_cast<int>(T/M_PI); 

    // T for 0 revs and moderate elliptical orbits
    double T_00 = acos(lambda) + lambda * sqrt(1 - lambda_squared);

    // multiple revs but moderate elliptical orbits
    double T_0 = (T_00 + M_max * M_PI); 
    // T at parabolic limit 
    double T_1 = 2.0/3.0 * (1.0 - lambda_cubed); 

    /*
    
    std::cout << "lambda cube : " << lambda_cubed << std::endl;
    std::cout << "T_1 : " << T_1 << std::endl;
    std::cout << "T_00 : " << T_00 << std::endl;
    std::cout << "T_0 : " << T_0 << std::endl;
   
    */

    // see equation 21 in D. Izzo "Revisiting Lambert's Problem"
    // initialize derivatives
    double first_derivative = 0.0; 
    double second_derivative = 0.0; 
    double third_derivative = 0.0;

    // if time of flight allows for multiple revolutions
    if (M_max > 0) {
        
        // if T < the orbital period for M_max moderate elliptical orbits  
        if (T < T_0){
        
            double T_min = T_0; 
            double x_now = 0.0; 
            double x_next = 0.0; 
            double error = 1.0;
            int iterations = 0;


            ///// START HALLEY ITERATIONS /////

            // ITERATION PROCEDURE from R. H. Gooding, "A Procedure for the Solution of Lambert's Orbital Boundary-Value Problem" (1970)

             while (1){
        
                taking_derivatives(first_derivative, second_derivative, third_derivative, lambda, x_now, T_min);
        
                if (first_derivative != 0.0) {

                    x_next = x_now - first_derivative * second_derivative/((second_derivative*second_derivative) - first_derivative * third_derivative / 2.0);

                }

                error = fabs(x_now - x_next);

                if (error < 1e-13 || iterations > 12){
                    break;
                }

                calc_time_of_flight(T_min, x_next, M_max, lambda);
                x_now = x_next;
                iterations++;

            }

            // if the minimum time of flight is > requested, then subtract a rev
            if (T_min > T){

                M_max = M_max - 1;
            }

        }
    
    }

    // minimum of possible and requested revolutions
    M_max = std::min(max_revs, M_max);

    /// Hereafter directly adopted from pykep ///
    // create output variables
    std::vector<double> m_v1;
    std::vector<double> m_v2;
    std::vector<int> m_iters;
    std::vector<double> m_x;

    // allocate memory  )
    m_v1.resize(static_cast<size_t>(M_max) * 2 + 1);
    m_v2.resize(static_cast<size_t>(M_max) * 2 + 1);
    m_iters.resize(static_cast<size_t>(M_max) * 2 + 1);
    m_x.resize(static_cast<size_t>(M_max) * 2 + 1);

    // 3 - We may now find all solutions in x,y
   
    // if the requested tof is longer than a less than single rev moderate elliptical orbit
    if (T >= T_00) {
        m_x[0] = -(T - T_00) / (T - T_00 + 4.0);
    } 
    // if not, then if requested tof is less than or equal to tof at parabolic limit 
    else if (T <= T_1) {
        m_x[0] = T_1 * (T_1 - T) / (2.0 / 5.0 * (1 - lambda_cubed * lambda_squared) * T) + 1.0;
    } 
    
    // if the tof requested is longer than that at the parabolic limit but less than that of a single rev moderate elliptical orbit
    else {
        m_x[0] = pow((T / T_00), 0.69314718055994529 / log(T_1 / T_00)) - 1.0;
    }

    //std::cout << "before householder: " << m_x[0] << " tof :" << T_00 <<" , " << T_1 << std::endl;


    // 3.1.2 Householder iterations
    m_iters[0] = iterate_householder(T, m_x[0], 0, 1e-5, 15,lambda);
    
    // 3.2 multi rev solutions
    double tmp;


    for (decltype(M_max) i = 1; i < M_max + 1; ++i) {
        // 3.2.1 left Householder iterations
        tmp = pow((i * M_PI + M_PI) / (8.0 * T), 2.0 / 3.0);
        m_x[2 * i - 1] = (tmp - 1) / (tmp + 1);
        m_iters[2 * i - 1] = iterate_householder(T, m_x[2 * i - 1], i, 1e-8, 15, lambda);
        // 3.2.1 right Householder iterations
        tmp = pow((8.0 * T) / (i * M_PI), 2.0 / 3.0);
        m_x[2 * i] = (tmp - 1) / (tmp + 1);
        m_iters[2 * i] = iterate_householder(T, m_x[2 * i], i, 1e-8, 15, lambda);
    }



    //std::cout << "xlist constructed with elements : " << m_x[0] << std::endl;



    double gamma = sqrt(mu * s / 2.0); 
    double rho = (r1_norm - r2_norm)/ c ; 
    double sigma = sqrt(1 - (rho*rho)); 
    double v_r1, v_t1, v_r2, v_t2, y, x_loop, y_loop;
    
    std::vector<arma::vec> velocity_solutions;

    //std::cout << "shape of dirvecs r1, r2 :" << dirvec_r1.n_elem << " , " << dirvec_r2.n_elem << std::endl;
    //std::cout << "shape of dirvecs t1, t2 :" << dirvec_t1.n_elem << " , " << dirvec_t2.n_elem << std::endl;

    // looping over all xy pairs returned from findxy
    for (int elem = 0; elem < m_x.size(); ++elem){

        x_loop = m_x[elem];
        double y_arg = 1.0 - (lambda_squared) + (lambda_squared) * (x_loop*x_loop);

        if (!(y_arg >= 0.0)) {
            std::cerr << "lambert_solver: y_arg is negative (" << y_arg << ") for x=" << x_loop << ", lambda=" << lambda << ", skipping solution\n";
            // push NaNs to keep return size consistent (caller should check for NaNs)
            arma::vec nanvec = arma::vec(6);
            nanvec.fill(arma::datum::nan);
            velocity_solutions.push_back(nanvec);
            continue;
        }

        y_loop = sqrt(y_arg);

        // diagnostics for NaN sources
        if (r1_norm == 0.0 || r2_norm == 0.0) {
            std::cerr << "lambert_solver: zero radius encountered r1_norm=" << r1_norm << " r2_norm=" << r2_norm << "\n";
        }
        double rho_check = (c == 0.0) ? arma::datum::nan : (r1_norm - r2_norm) / c;
        double sigma_check = (rho_check != rho_check) ? arma::datum::nan : sqrt(1 - (rho_check*rho_check));
        if (!std::isfinite(rho_check) || !std::isfinite(sigma_check)) {
            std::cerr << "lambert_solver: invalid rho/sigma rho_check=" << rho_check << " sigma_check=" << sigma_check << "\n";
        }

        v_r1 = gamma*( (lambda*y_loop - x_loop) -  rho*( lambda*y_loop + x_loop) )  / r1_norm;
        v_r2 = -gamma*( (lambda*y_loop - x_loop) +  rho*( lambda*y_loop + x_loop) )  / r2_norm;

        v_t1 = gamma*sigma*( y_loop + lambda*x_loop ) / r1_norm;
        v_t2 = gamma*sigma*(y_loop + lambda*x_loop) / r2_norm;

        arma::vec v1 = v_r1*dirvec_r1 + v_t1*dirvec_t1;
        arma::vec v2 = v_r2*dirvec_r2 + v_t2*dirvec_t2;

        velocity_solutions.push_back(arma::join_cols(v1,v2));

    }

    return velocity_solutions;


}


arma::vec3 get_v1(std::vector<arma::vec> solutions, int const solution_number){

    return solutions[solution_number].subvec(0,2);

}

arma::vec3 get_v2(std::vector<arma::vec> solutions, int const solution_number){

    return solutions[solution_number].subvec(3,5);

}

