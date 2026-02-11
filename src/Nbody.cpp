#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <armadillo>
#include <chrono>
#include <random>

#include "Nbody.hpp"
#include "Local_search.hpp"
#include "data_access_lib.hpp"

/*
Author : S. Nigudkar (2025)

Functions that allow the initialization and propagation of a satellite constellation and a
service station in circular orbit. Available constellation arrangements are currently limited
to - Delta Walker.

Units used - m,s,kg

Time Integrator used - Runge-Kutta 4

*/

// sim functions //

double deg_to_rads(double deg){

    //converts measured angles from degrees to radians

    return deg * M_PI / 180.0;

}



double get_inclination(arma::vec r, arma::vec v){

    double r_norm = arma::norm(r);
    double v_norm = arma::norm(v);

    // angular momentum
    arma::vec h = arma::cross(r,v);
    double h_norm = arma::norm(h);

    // inclination
    double i = std::acos(h(2)/h_norm);

    return i;
}





orbital_elements orb_elems_from_rv(arma::vec r, arma::vec v){


    // a, e, i, RAAN, true anomaly and augment_of_periapsis
    orbital_elements orb_elems;

    double r_norm = arma::norm(r);
    double v_norm = arma::norm(v);

    // angular momentum
    arma::vec h = arma::cross(r,v);
    double h_norm = arma::norm(h);

    // inclination
    double i = std::acos(h(2)/h_norm);

    // RAAN
    arma::vec k = {0,0,1};

    arma::vec n = arma::cross(k,h);

    double n_norm = arma::norm(n);


    double RAAN;

    if (n_norm < 1e-10) {

      RAAN = 0.0;  // Convention
    }

    else if (n(1) >= 0 ){

        RAAN =  std::acos(n(0)/n_norm);

    }

    else{

        RAAN = 2 * M_PI - std::acos(n(0)/n_norm);

    }


    // eccentricity


    arma::vec e = ((v_norm*v_norm - MU_EARTH/r_norm) * r - arma::dot(r,v)*v) / MU_EARTH;


    double e_norm = arma::norm(e);


    // augment_of_periapsis

    double p_arg;


    if (e_norm < 1e-10 || n_norm < 1e-10) {
      p_arg = 0.0;  // Convention
    }



    else if (e(2) >= 0){

        p_arg = std::acos(arma::dot(n,e)/ (n_norm*e_norm));

    }


    else{

        p_arg = 2*M_PI - std::acos(arma::dot(n,e)/ (n_norm*e_norm));

    }



    // True anomaly

    double nu;

    if (e_norm < 1e-10) {
      nu = 0.0;  // Convention
    }

    else if (arma::dot(r,v) >= 0){

      nu = std::acos(arma::dot(e,r) / (e_norm*r_norm));

    }

    else{

      nu = 2*M_PI - std::acos(arma::dot(e,r) / (e_norm*r_norm));

    }




    // semi major axis

    double energy =  (0.5 * v_norm*v_norm) - MU_EARTH/r_norm;
    double a;


    // when the energy approaches 0 we have a parabolic orbit
    if (std::abs(energy) < 1e-10){

        //for parabolic orbits energy = 0 so 'a' tends to infinity

        a = std::numeric_limits<double>::infinity();

    }

    else{

        // for elliptical orbits e<0 and for hyperbolic orbits e>0
        // so a can be finite

        a = -MU_EARTH/ (2.0 * energy);


    }



    // we have now calculated all orbital elements
    // now assign these to the orbital_elements object


    orb_elems.inclination = i;
    orb_elems.RAAN = RAAN;
    orb_elems.augment_of_periapsis = p_arg;
    orb_elems.true_anomaly = nu;
    orb_elems.semi_major_axis = a;
    orb_elems.eccentricity = e_norm;


    return orb_elems;

}








arma::vec3 acceleration_due_to_central_body (const arma::vec3& r){

    double r_magnitude = arma::norm(r);

    if(r_magnitude == 0){

        return arma::vec3();
    }

    return r * (-MU_EARTH / (pow(r_magnitude,3)));
}

// acceleration taking j2 perturbations into account
arma::vec3 acceleration_due_to_J2(const arma::vec3& r){

    double rx = r(0);
    double ry = r(1);
    double rz = r(2);

    double r2 = rx*rx + ry*ry + rz*rz;
    double r5 = pow(r2, 2.5);

    double z2 = rz*rz;

    double factor = 1.5 * J2 * MU_EARTH * R_EARTH*R_EARTH / r5;

    double common = 5.0 * z2 / r2 - 1.0;

    arma::vec3 a;

    a(0) = rx * factor * common;
    a(1)= ry * factor * common;
    a(2) = rz * factor * (5.0 * z2 / r2 - 3.0);
    return a;
}


// compute accelerations given froce options
arma::vec3 compute_acceleration(const satellite_object& sat, const std::vector<satellite_object>& all_sats, const force_model& force_options){

    arma::vec3 acceleration_total = acceleration_due_to_central_body(sat.r);

    // adding a J2 perturbation
    // (Non-keplerian)

    if(force_options.includeJ2){
        acceleration_total += acceleration_due_to_J2(sat.r);
    }


    // the perturbations in orbits caused by the mutual attraction between satellites
    // (Non keplerian)
    if(force_options.includeMutual){

        // calculate the mutual force applied by all sats
        for(const auto& other_satellite : all_sats){
            // exclude self from calculation (force applied by satellite on itself)
            if(&other_satellite == &sat) continue;

            // distance between satellites
            arma::vec3 dr_sat = other_satellite.r - sat.r;

            double d = arma::norm(dr_sat);

            if(d == 0) continue;

                acceleration_total += dr_sat * (G_CONST * other_satellite.mass / std::pow(d,3));
        }
    }
    return acceleration_total;
}


// 4th order runge kutta integration to propogate satellite objects forwards in time.
void runge_kutta_step(std::vector<satellite_object>& sats, double dt, const force_model& force_options){

    int num_sats = sats.size();

    // initialize vectors to store RK additions for both position and velocity
    std::vector<arma::vec> k1r(num_sats), k1v(num_sats), k2r(num_sats), k2v(num_sats), k3r(num_sats), k3v(num_sats), k4r(num_sats), k4v(num_sats);

    // we will store the outputs of the rk4 passes in this vector
    std::vector<satellite_object> temp = sats;

    // iterating through all satellites

    // RK4 Algo

    // first pass = f(yn)
    // k1 = f(yn + h * (fist pass) / 2 )
    // k2 = f(yn + h * (k1) / 2 )
    // k3 = f(yn + h * k2 )
    // k4 = y_(n+1) = y_n + h/6 (first pass + 2k1 + 2k2 + k3)

    // first pass is already stored in sats

    // now we use the values in sats to calculate the rk1
    for(int i=0;i<num_sats;i++) k1r[i] = sats[i].v;
    for(int i=0;i<num_sats;i++) k1v[i] = compute_acceleration(sats[i], sats, force_options);
    for(int i=0;i<num_sats;i++){

        temp[i].r = sats[i].r + k1r[i]*(0.5*dt);
        temp[i].v = sats[i].v + k1v[i]*(0.5*dt);
    }

    // rk2

    for(int i=0;i<num_sats;i++) k2r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k2v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){

        temp[i].r = sats[i].r + k2r[i]*(0.5*dt);
        temp[i].v = sats[i].v + k2v[i]*(0.5*dt);
    }


    // rk3

    for(int i=0;i<num_sats;i++) k3r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k3v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){

        temp[i].r = sats[i].r + k3r[i]*dt;
        temp[i].v = sats[i].v + k3v[i]*dt;
    }


    // using above passes to compute y_n+1

    for(int i=0;i<num_sats;i++) k4r[i] = temp[i].v;
    for(int i=0;i<num_sats;i++) k4v[i] = compute_acceleration(temp[i], temp, force_options);
    for(int i=0;i<num_sats;i++){

        sats[i].r = sats[i].r + (k1r[i] + k2r[i]*2.0 + k3r[i]*2.0 + k4r[i]) * (dt/6.0);
        sats[i].v = sats[i].v + (k1v[i] + k2v[i]*2.0 + k3v[i]*2.0 + k4v[i]) * (dt/6.0);
    }



}


// distributing satellites in a walker constellation
std::vector<satellite_object> build_walker_constellation(int num_planes, int total_satnum, int phase ,
                                             double altitude_m,
                                             double inclination_rad,
                                             double sat_mass) {

    std::vector<satellite_object> sats;

    // predecided size of array
    sats.reserve(total_satnum);


    double a = altitude_m;
    // num satellites per plane
    int S = total_satnum / num_planes ;
    // the right accession (a total of 360 deg or 2pi rads) is divided up into n orbital planes.
    double RAAN_step = 2.0 * M_PI / num_planes;

    // for each plane,
    for (int p = 0; p < num_planes; ++p) {

        double RAAN = p * RAAN_step;

        // for each sat in this plane
        for (int s = 0; s < S; ++s) {

            // apply phasing
            double phase = 2.0 * M_PI * ( (s + (double)phase * p / num_planes) / S );

            orbital_elements el;
            el.semi_major_axis = a;
            el.eccentricity = 0.0;
            el.inclination = inclination_rad;
            el.RAAN = RAAN;
            el.augment_of_periapsis= 0.0;
            el.true_anomaly = phase;
            std::string name = "sat_" + std::to_string(p*S + s);
            sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords(name, el, sat_mass));
        }
    }
    return sats;
}





std::vector<satellite_object> circular_orbits(int total_satnum,
                                             double altitude_m,
                                             double inclination_rad,
                                             double RAAN , double satmass){


    std::vector<satellite_object> sats;
    // predecided size of array
    sats.reserve(total_satnum);

    double a = altitude_m;

    orbital_elements el;


    el.semi_major_axis = a;
    el.eccentricity = 0.0;
    el.inclination = inclination_rad;
    el.RAAN = RAAN;

    double phase = 2 * M_PI / total_satnum ;

    el.augment_of_periapsis= 0.0;

    for (int p = 0; p < total_satnum; p++){

        el.true_anomaly = phase * p;

        std::string name = "sat_"+std::to_string(p) ;

        sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords(name, el, satmass));

    }


    return sats;

}




void showProgressBar(int progress, int total, int barWidth ) {
    float ratio = static_cast<float>(progress) / total;
    int pos = static_cast<int>(barWidth * ratio);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "█";   // filled part
        else std::cout << "░";           // empty part
    }
    std::cout << "] " << std::setw(3) << int(ratio * 100) << "%\r";
    std::cout.flush();
    std::cout<<"\n";
}





void run_simulation(std::string save_to_file, double t_final, force_model fmodel){

    std::vector<satellite_object> sats;

    sats = build_walker_constellation(num_planes, num_clients, phase,
                                      a_client,inclination,1.0);

    double period = 2*M_PI * sqrt(pow(a_client,3) / MU_EARTH);

    //init service satellite in circular orbit

    orbital_elements svc;
    svc.semi_major_axis = a_depot;
    svc.eccentricity = 0.001;
    svc.inclination = inclination;
    svc.RAAN = 0.0;
    svc.augment_of_periapsis = 10.0 * M_PI/180.0;
    svc.true_anomaly = 0.0;

    sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords("service_1", svc, 1.0));

    //std::cout << "# time(s) index name x y z vx vy vz\n";

    //write header
    std::ofstream csv("../data/"+save_to_file);
    csv << "time_s,index,name,x,y,z,vx,vy,vz,semi_major_axis,eccentricity,inclination,RAAN,arg_periapsis,true_anomaly\n";


    double t = 0.0;
    int step = 0;


    std::cout<<"Running Simulation";
    while(t < (t_final - 1e-9)){
        // only write outputs to csv files every 10 timesteps
        // I am just thinning the datafile here so we can store it easily

        //if(step % 10 == 0){

            for(size_t i=0;i<sats.size();i++){

                // Calculate orbital elements for current state
                orbital_elements elems = orb_elems_from_rv(sats[i].r, sats[i].v);

                csv << t << "," << i << "," << sats[i].satname << ","
                    << sats[i].r(0) << "," << sats[i].r(1) << "," << sats[i].r(2) << ","
                    << sats[i].v(0) << "," << sats[i].v(1) << "," << sats[i].v(2) << ","
                    << elems.semi_major_axis << "," << elems.eccentricity << ","
                    << elems.inclination << "," << elems.RAAN << ","
                    << elems.augment_of_periapsis << "," << elems.true_anomaly << "\n";

            }

            //}

        runge_kutta_step(sats, tstep_size, fmodel);

        t += tstep_size;
        step++;
        //showProgressBar(t, t_final);


    }

    std::cout << std::endl << "Done!" << std::endl;

    csv.close();
    std::cout << "CSV saved to ../data/"+save_to_file+"\n";


}



void run_simulation(    std::string save_to_file, std::string arrangement, double t_final,double dt,
                        double altitude_m, double num_planes, double num_satellites,
                        double relative_phase, double inclination_in_radians,force_model fmodel,
                        double satmass,double idiff){


    double altitude = R_EARTH + altitude_m;

    double period = 2*M_PI * sqrt(pow(altitude,3) / MU_EARTH);

    std::cout << "orbital period: " <<  period;

    std::vector<satellite_object> sats;

    if (arrangement == "walker_delta"){

        sats = build_walker_constellation(num_planes, num_satellites, relative_phase,
                                                                        altitude,
                                                                        56.0 * M_PI/180.0,
                                                                        satmass);

    }

    if (arrangement == "circular_orbits"){


        sats = circular_orbits(num_satellites,
                                altitude,inclination_in_radians,
                                0.0, satmass);



    }

    if (arrangement == "flower_constellation"){




        //std::vector<satellite_object> sats = build_coplanar_circular_orbit();

    }

    if (arrangement == "spiral_constellation"){

        //std::vector<satellite_object> sats = build_coplanar_circular_orbit();


    }

    if (arrangement == "custom"){

        //std::vector<satellite_object> sats = build_coplanar_circular_orbit();

    }


    //init service satellite in circular orbit

    orbital_elements svc;
    svc.semi_major_axis = R_EARTH + 720000.0;
    svc.eccentricity = 0.001;
    svc.inclination = inclination_in_radians + idiff; //48.0 * M_PI/180.0;
    svc.RAAN = 0.0;  //30.0 * M_PI/180.0;
    svc.augment_of_periapsis = 10.0 * M_PI/180.0;
    svc.true_anomaly = 0.0;
    sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords("service_1", svc, 1.0));

    //std::cout << "# time(s) index name x y z vx vy vz\n";

    //write header
    std::ofstream csv("../data/"+save_to_file);
    csv << "time_s,index,name,x,y,z,vx,vy,vz,semi_major_axis,eccentricity,inclination,RAAN,arg_periapsis,true_anomaly\n";


    double t = 0.0;
    int step = 0;


    std::cout<<"Running Simulation";
    while(t < (t_final - 1e-9)){
        // only write outputs to csv files every 10 timesteps
        // I am just thinning the datafile here so we can store it easily

        if(step % 10 == 0){

            for(size_t i=0;i<sats.size();i++){

                // Calculate orbital elements for current state
                orbital_elements elems = orb_elems_from_rv(sats[i].r, sats[i].v);

                csv << t << "," << i << "," << sats[i].satname << ","
                    << sats[i].r(0) << "," << sats[i].r(1) << "," << sats[i].r(2) << ","
                    << sats[i].v(0) << "," << sats[i].v(1) << "," << sats[i].v(2) << ","
                    << elems.semi_major_axis << "," << elems.eccentricity << ","
                    << elems.inclination << "," << elems.RAAN << ","
                    << elems.augment_of_periapsis << "," << elems.true_anomaly << "\n";

            }

        }

        runge_kutta_step(sats, dt, fmodel);

        t += dt;
        step++;
        //showProgressBar(t, t_final);


    }

    std::cout << std::endl << "Done!" << std::endl;

    csv.close();
    std::cout << "CSV saved to ../data/WalkerDelta.csv\n";



}



bool is_in_vec(const std::vector<double>& v, double x) {
    return std::find(v.begin(), v.end(), x) != v.end();
}




struct constellation_struct{
                        int num_planes;
                        int num_satellites;
                        int relative_phase;
                        double altitude;
                        double inclination;
                        double satmass;

};



class Simulation{


    public:

        force_model force_options;
        int numsats;
        std::vector<satellite_object> sats;
        std::vector<double> taken_inclinations;
        std::vector<double> taken_a;
        // constructor
        Simulation() : force_options(true, false) {
            numsats = 0;

        };


        Simulation(force_model force_opts) : force_options(force_opts) {
            numsats=0;

        };



        void add_sat(orbital_elements el, std::string name, double satmass){


                sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords(name, el, satmass));


        }


        void add_sats(orbital_elements el,int totalsats,double phase_lim,double satmass){


            if(is_in_vec(taken_a,el.semi_major_axis) && is_in_vec(taken_inclinations, el.inclination)){

                std::cout<<"The desired orbit is already populated, please choose different (a,i)"<<std::endl;

            }

            else{


                for (int i= 0 ; i < totalsats ; i++){

                    numsats++;

                    orbital_elements el_in = el;

                    double phase = 2 * M_PI/phase_lim *  i / totalsats ;

                    el_in.true_anomaly = phase;

                    sats.push_back(satellite_object::from_satellite_normal_to_ECI_coords("sat_"+std::to_string(numsats), el_in, satmass));

                }


                taken_inclinations.push_back(el.inclination);
                taken_a.push_back(el.semi_major_axis);

            }

        };



        void add_constellation(std::string constellation_type,constellation_struct constellation){


            if (numsats > 0){

                std::cout << "unable to populate constellation because orbits already exist"<<std::endl;
                std::cout << "please start with an empty simulation environment and initialize constellation before other orbits"<<std::endl;

            }

            else{


            if (constellation_type == "walker_delta"){

                std::vector<satellite_object> constellation_array = build_walker_constellation( constellation.num_planes,
                                                                                                constellation.num_satellites,
                                                                                                constellation.relative_phase,
                                                                                                constellation.altitude,
                                                                                                constellation.inclination,
                                                                                                constellation.satmass);



                sats = constellation_array;

                taken_inclinations.push_back(constellation.inclination);
                taken_a.push_back(constellation.altitude);

                numsats+=constellation.num_satellites;


            }

            }


        };


        // destructor
        ~Simulation();

        // now we define functions to actually run the simulation

        void run_simulation(double t_end,double dt, std::string file_to_write){

            std::ofstream csv("../data/"+file_to_write);
            csv << "time_s,index,name,x,y,z,vx,vy,vz,semi_major_axis,eccentricity,inclination,RAAN,arg_periapsis,true_anomaly\n";

            double t = 0.0;
            std::cout<<"Running Simulation";
            int step=0;

            while(t < (t_end - 1e-9)){
                // only write outputs to csv files every 10 timesteps
                // I am just thinning the datafile here so we can store it easily

                if(step % 10 == 0){

                    for(size_t i=0;i<sats.size();i++){

                        // Calculate orbital elements for current state
                        orbital_elements elems = orb_elems_from_rv(sats[i].r, sats[i].v);

                        csv << t << "," << i << "," << sats[i].satname << ","
                            << sats[i].r(0) << "," << sats[i].r(1) << "," << sats[i].r(2) << ","
                            << sats[i].v(0) << "," << sats[i].v(1) << "," << sats[i].v(2) << ","
                            << elems.semi_major_axis << "," << elems.eccentricity << ","
                            << elems.inclination << "," << elems.RAAN << ","
                            << elems.augment_of_periapsis << "," << elems.true_anomaly << "\n";

                    }

                }

                runge_kutta_step(sats, dt, force_options);

                t += dt;
                step++;
                //showProgressBar(t, t_final);


            }

            std::cout << std::endl << "Done!" << std::endl;

            csv.close();
            std::cout << "CSV saved to ../data/WalkerDelta.csv\n";


        };


    private:

        bool check_for_collisions(orbital_elements el){

            return false;

        }




};



schedule_struct create_instance(int num_visits, DataFrame simfile){

    std::vector<double> t_depart;
    std::vector<double> t_arrive;

    double initdeltav = 1e22;

    double service_time = 18000;

    double t_final = 350000 * num_visits;

    init_dep_arrival_times_strict_timespan( t_depart, t_arrive, t_final, service_time, num_visits );

    SatelliteData sat_data(simfile);

    std::vector<std::string> satnames;
    std::vector<std::string> service_names;
    std::vector<std::string> client_names;

    for(std::string name : sat_data.names ){

        if(name.find("service") != std::string::npos){
            service_names.push_back(name);
        }

        else{

            client_names.push_back(name);

        }

    }

    std::mt19937 gen(std::random_device{}());
    std::bernoulli_distribution service_dist(0.3);
    std::bernoulli_distribution client_dist(0.7);
    std::uniform_int_distribution<size_t> service_id_dist(0, service_names.size() - 1);


    int sats = 0 ;


    while(sats < num_visits){

        //bool value_service = service_dist(gen);
        //bool value_client  = client_dist(gen);

        //if(value_service == true){

        size_t idx = service_id_dist(gen);
        satnames.push_back(service_names[idx]);


            //}

        //if(value_client == false){

            std::uniform_int_distribution<size_t> client_id_dist(0, service_names.size() - 1);
            idx = client_id_dist(gen);

            satnames.push_back(client_names[idx]);
            client_names.erase(client_names.begin() + idx);

            sats++;
        //}

    }

    size_t idx = service_id_dist(gen);
    satnames.push_back(service_names[idx]);

    schedule_struct schedule_base = create_schedule_lambert_only(initdeltav, t_arrive, t_depart, satnames, simfile, service_time);
    return schedule_base;


}
