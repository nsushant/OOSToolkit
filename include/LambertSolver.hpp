

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



void taking_derivatives(    double &DT,  double &DDT, double &DDDT, 
                            double const lambda, const double x, 
                            const double T  );





void tof_lagrange(double &tof, const double x, const int N, const double lambda);


double hypergeometricF(double z, double tol);


void calc_time_of_flight(   double &time_of_flight,const double x, 
                            const int N , const double lambda    );



int iterate_householder(const double T, double &x0, 
                        const int N, const double eps, 
                        const int iter_max,const double lambda);




// the lambert solver is going to return velocitiies v1, v2 which have to be applied at t1,t2 respectively 
// inputs - r1 (t1) , r2 (t2), time of flight (t2-t1), mu, direction (1 for prograde, -1 for retrograde), max_revolutions (default 0)
std::vector<arma::vec> lambert_solver(  const arma::vec3 r1, const arma::vec3 r2,const 
                                        double time_of_flight, const double mu ,
                                        const int retrograde, const int max_revolutions);



arma::vec3 get_v1(std::vector<arma::vec> solutions, int const solution_number);

arma::vec3 get_v2(std::vector<arma::vec> solutions, int const solution_number);

