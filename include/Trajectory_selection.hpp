
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <armadillo>
#include "LambertSolver.hpp"
#include "data_access_lib.hpp"
#include "Nbody.hpp"
#include "Local_search.hpp"

// Forward declarations for DP types
struct lookupkey;
struct HashKey;


/*
Author : S. Nigudkar (2025)

The aim of the below function is to choose lambert parameters r1,r2,tof and compute the associated delta Vs
for the most efficient orbital transfers given a requested service deadline t_service, the start time
t_prev, the name of the client satellite (departure destination) and that of the service satellite (arrival destination).

*/

// Helper struct to consolidate DataFrame access patterns
struct SatelliteData {
    arma::vec x, y, z, vx, vy, vz, t;
    std::vector<std::string> names;

    explicit SatelliteData(DataFrame& simfile) {
        x = simfile.getNumeric("x");
        y = simfile.getNumeric("y");
        z = simfile.getNumeric("z");
        vx = simfile.getNumeric("vx");
        vy = simfile.getNumeric("vy");
        vz = simfile.getNumeric("vz");
        t = simfile.getNumeric("time_s");
        names = simfile["name"];
    }
};

inline arma::uvec find_satellite_indices(const SatelliteData& data, const std::string& satname) {
    return find_idxs_of_match(data.names, satname);
}

//class to get orbital elems.

class OrbitalElementAccessor {
private:
    DataFrame* simfile;

    // Cache for frequently accessed (sat_name, time) combinations
    static std::unordered_map<std::string, std::map<double, orbital_elements>> elem_cache;

public:
    OrbitalElementAccessor(DataFrame* df) : simfile(df) {}

    orbital_elements get_elements_at_time(const std::string& sat_name, double time) {
        // Check local cache first
        auto sat_it = elem_cache.find(sat_name);
        if (sat_it != elem_cache.end()) {
            auto time_it = sat_it->second.find(time);
            if (time_it != sat_it->second.end()) {
                return time_it->second;
            }
        }

        // Load from DataFrame (same pattern as current position/velocity access)
        std::vector<std::string> satnames = (*simfile)["name"];
        arma::uvec sat_idxs = find_idxs_of_match(satnames, sat_name);

        arma::vec t_sat = simfile->getNumeric("time_s");
        arma::uvec time_idx = arma::find(t_sat == time);

        // Extract orbital elements from enhanced CSV
        arma::vec a_sat = simfile->getNumeric("semi_major_axis");
        arma::vec e_sat = simfile->getNumeric("eccentricity");
        arma::vec i_sat = simfile->getNumeric("inclination");
        arma::vec raan_sat = simfile->getNumeric("RAAN");
        arma::vec omega_sat = simfile->getNumeric("arg_periapsis");
        arma::vec nu_sat = simfile->getNumeric("true_anomaly");

        orbital_elements elems;
        elems.semi_major_axis = a_sat(time_idx(0));
        elems.eccentricity = e_sat(time_idx(0));
        elems.inclination = i_sat(time_idx(0));
        elems.RAAN = raan_sat(time_idx(0));
        elems.augment_of_periapsis = omega_sat(time_idx(0));
        elems.true_anomaly = nu_sat(time_idx(0));

        // Cache result
        elem_cache[sat_name][time] = elems;

        return elems;
    }
};



void find_optimal_trajectory(std::string service_satname, std::string client_satname, double &t_prev,
    double &t_request, arma::vec &v1sol, arma::vec &v2sol, arma::vec &r1sol, arma::vec &r2sol,
    double &tof_optimal, std::vector<arma::vec> &trajs, DataFrame simfile, std::string method,
    bool write_to_file, double &DeltaVMinima);

void find_optimal_trajectory(std::string service_satname, std::string client_satname, double &t_prev,
    double &t_request, arma::vec &v1sol, arma::vec &v2sol, arma::vec &r1sol, arma::vec &r2sol,
    double &tof_optimal, std::vector<arma::vec> &trajs, DataFrame simfile, std::string method,
    bool write_to_file, double &DeltaVMinima, std::unordered_map<lookupkey, double, HashKey>& lookup);

double compact_optimal_calc(std::string satname1, std::string satname2, double departure_time, double arrival_time, DataFrame simfile);

void find_optimal_trajectory_no_iter(std::string service_satname, std::string client_satname, double t_departure,
    double t_arrival, DataFrame simfile, double &DeltaVMinima, std::string method="edelbaum");

void run_exhaustive_search( std::string sat_from, std::string sat_to, double &t_from, double &t_to, double &deltaV_change, std::string simfilename,
                            std::string outputfilename, std::string method = "edelbaum" );
