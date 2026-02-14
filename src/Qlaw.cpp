#include <armadillo>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <iostream>
#include <algorithms>
#include <vector>
#include <cmath>
#include <cmath>
#include "Nbody.hpp"


//equinoctual orbital elems
double calcp(double a,double e){

    return a*(1-e*e);
}

double calcf(double e, double argperi, double raan){

    return e*std::cos(argperi + raan);
}

double calcg(double e, double argperi, double raan){

    return e*std::sin(argperi + raan);
}

double calch(double i, double raan){

    return std::tan(i/2)*std::cos(raan);
}


double calck(double i, double raan){

    return std::tan(i/2)*std::sin(raan);
}


double calcl(double raan, double argperi, double true_anomaly){

    return argperi + raan + true_anomaly;

}


struct EQorbital_elements{

    double a,p,f,g,h,k,L;

};

struct optparams{

    double wa,wf,wg,wh,wk,wl,wscl;

    double eta_a,eta_r,m,n,r,rot_angle;

};


EQorbital_elements convert_kep_to_eq(orbital_elements oe){

    EQorbital_elements eq_oe;

    eq_oe.p = calcp(oe.semi_major_axis, oe.eccentricity);

    eq_oe.f = calcf(oe.eccentricity, oe.augment_of_periapsis, oe.RAAN);

    eq_oe.g = calcg(oe.eccentricity, oe.augment_of_periapsis, oe.RAAN);

    eq_oe.h = calch(oe.inclination, oe.RAAN);

    eq_oe.k = calck(oe.inclination, oe.RAAN);

    eq_oe.L = calcl(oe.RAAN,  oe.augment_of_periapsis, oe.true_anomaly);

    eq_oe.a = oe.semi_major_axis;

    return eq_oe;

}




double aug_oe(EQorbital_elements eq_oeC, EQorbital_elements eq_oeT, optparams weights, double rpmin){

    double deltaL= eq_oeC.L - eq_oeT.L;

    return ((2*weights.wl)/M_PI) *( eq_oeT.a - ( rpmin )/(1 - std::sqrt((eq_oeC.f*eq_oeC.f) +
           (eq_oeC.g*eq_oeC.g))))*std::atan( weights.wscl*deltaL) + eq_oeT.a;


}



struct derivatives {

    double amax,fmax,gmax,hmax,kmax;

};



derivatives calcdt(double F_accel, EQorbital_elements current_oe){

    double s = 1 + current_oe.k * current_oe.k + current_oe.h* current_oe.h;

    derivatives ddt;

    ddt.fmax = 2*F_accel*std::sqrt(current_oe.p/MU_EARTH);
    ddt.gmax = 2*F_accel*std::sqrt(current_oe.p/MU_EARTH);
    ddt.hmax = 0.5 * F_accel * std::sqrt(current_oe.p/MU_EARTH) * s*s / (std::sqrt(1-current_oe.g*current_oe.g) + current_oe.f);
    ddt.kmax = 0.5 * F_accel * std::sqrt(current_oe.p/MU_EARTH) * s*s / (std::sqrt(1-current_oe.f*current_oe.f) + current_oe.g);
    ddt.amax = 2*F_accel*current_oe.a*std::sqrt(current_oe.a/MU_EARTH)*std::sqrt((1 + std::sqrt(current_oe.f*current_oe.f + current_oe.g*current_oe.g))/(1 - std::sqrt(current_oe.f*current_oe.f + current_oe.g*current_oe.g)) );

    return ddt;

}




//Varga and Perez Q-law
double Qval(orbital_elements oec, orbital_elements oet, double r_pmin,
            double k, double wp, optparams weights, double F_accel){

    EQorbital_elements current_oe = convert_kep_to_eq(oec);

    EQorbital_elements target_oe = convert_kep_to_eq(oet);

    //target_oe.a = aug_oe(current_oe, target_oe, weights, r_pmin);

    double r_p = current_oe.p/(1+oec.eccentricity);

    double P = std::exp(k*(1-r_p/r_pmin));

    double penalty = (1 + wp*P);

    double s_a = std::pow((1+ std::pow((std::abs(current_oe.a - target_oe.a) / (weights.m*target_oe.a)), weights.n)),1/weights.r);

    double summation = 0.0;

    auto sumcalc = [](double s, double w, double oe_ddx, double oe1, double oe2){

        return s*w*(  (oe1-oe2) / oe_ddx  ) * (  (oe1-oe2) / oe_ddx  );
    };

    double s = 1 + current_oe.k * current_oe.k + current_oe.h* current_oe.h;

    derivatives ddt = calcdt(F_accel, current_oe);

    summation+= sumcalc(s_a, weights.wa, ddt.amax, current_oe.a, target_oe.a);
    summation+= sumcalc(1.0, weights.wf, ddt.fmax, current_oe.f, target_oe.f);
    summation+= sumcalc(1.0, weights.wg, ddt.gmax, current_oe.g, target_oe.g);
    summation+= sumcalc(1.0, weights.wh, ddt.hmax, current_oe.h, target_oe.h);
    summation+= sumcalc(1.0, weights.wk, ddt.kmax, current_oe.k, target_oe.k);

    return penalty*summation;

}



double dQdt(EQorbital_elements oe, derivatives ddt){





}
