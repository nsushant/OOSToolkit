#include <casadi/casadi.hpp>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include "Nbody.hpp"  // Include your header for tstep_size and t_final
#include "LowThrustAnalytical.hpp"

using namespace casadi;

// Structure to hold data from CSV
struct SatelliteState {
    double time;
    double x, y, z;
    double vx, vy, vz;
    orbital_elements orb;
};

class EdelbaumCallback : public Callback {
private:
    // Map: sat_index -> vector of states over time
    std::map<int, std::vector<SatelliteState>> sat_data;
    double min_time;
    double max_time;

public:
    EdelbaumCallback(const std::string& csv_file)
        : min_time(1e10), max_time(-1e10) {
        load_simulation_data(csv_file);
    }

    casadi_int get_n_in() override { return 4; }  // [sat_dep_idx, sat_arr_idx, t_dep, t_arr]
    casadi_int get_n_out() override { return 1; }  // [delta_v]

    std::vector<DM> eval(const std::vector<DM>& arg) const override {
        int sat_dep_idx = static_cast<int>(std::round(arg[0].scalar()));
        int sat_arr_idx = static_cast<int>(std::round(arg[1].scalar()));
        double t_dep = arg[2].scalar();
        double t_arr = arg[3].scalar();

        double transfer_time = t_arr - t_dep;

        // Enforce minimum time resolution using tstep_size from Nbody.hpp
        if (transfer_time < tstep_size) {
            return {1e10};  // Penalty
        }

        // Check time bounds (0 to t_final from Nbody.hpp)
        if (t_arr < 0 || t_arr > t_final || t_dep < 0 || t_dep > t_final) {
            return {1e10};  // Out of simulation range
        }

        // Check if satellites exist
        if (sat_data.find(sat_dep_idx) == sat_data.end() ||
            sat_data.find(sat_arr_idx) == sat_data.end()) {
            return {1e10};  // Invalid satellite
        }

        // Get orbital elements at ARRIVAL TIME for both satellites
        orbital_elements orb_departure = get_orbital_elements_at_time(sat_dep_idx, t_arr);
        orbital_elements orb_arrival = get_orbital_elements_at_time(sat_arr_idx, t_arr);

        // Call your Edelbaum function
        double dv = calculate_edelbaum_deltaV_orbital_elems(
            orb_departure,
            orb_arrival,
            transfer_time
        );

        return {dv};
    }

    // Tell CasADi we'll provide custom derivatives
    bool has_reverse(casadi_int nadj) const override { return nadj == 1; }

    // Implement reverse mode AD (gradient computation)
    std::vector<std::vector<DM>> call_reverse(const std::vector<DM>& arg,
                                               const std::vector<DM>& res,
                                               const std::vector<std::vector<DM>>& adj,
                                               bool always_inline) const {
        // Use finite differences for gradients
        double h = tstep_size;

        // Current function value
        double f0 = res[0].scalar();
        double adj_f = adj[0][0].scalar();  // Adjoint seed

        // Initialize gradients
        std::vector<DM> grad(4);

        // Gradient w.r.t. sat_dep (discrete - set to 0)
        grad[0] = 0.0;

        // Gradient w.r.t. sat_arr (discrete - set to 0)
        grad[1] = 0.0;

        // Gradient w.r.t. t_departure
        std::vector<DM> arg_dep = arg;
        arg_dep[2] = arg[2].scalar() + h;
        double f_dep = eval(arg_dep)[0].scalar();
        grad[2] = adj_f * (f_dep - f0) / h;

        // Gradient w.r.t. t_arrival
        std::vector<DM> arg_arr = arg;
        arg_arr[3] = arg[3].scalar() + h;
        double f_arr = eval(arg_arr)[0].scalar();
        grad[3] = adj_f * (f_arr - f0) / h;

        return {grad};
    }

    // Getter for time bounds
    double get_min_time() const { return min_time; }
    double get_max_time() const { return max_time; }

    // Get list of available satellite indices
    std::vector<int> get_satellite_indices() const {
        std::vector<int> indices;
        for (const auto& pair : sat_data) {
            indices.push_back(pair.first);
        }
        return indices;
    }

private:
    void load_simulation_data(const std::string& csv_file) {
        std::ifstream file(csv_file);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open CSV file: " + csv_file);
        }

        std::string line;
        std::getline(file, line);  // Skip header

        while (std::getline(file, line)) {
            std::istringstream ss(line);
            std::string token;

            SatelliteState state;
            int index;
            std::string name;

            // Parse CSV line: time_s,index,name,x,y,z,vx,vy,vz,semi_major_axis,
            //                 eccentricity,inclination,RAAN,arg_periapsis,true_anomaly

            std::getline(ss, token, ','); state.time = std::stod(token);
            std::getline(ss, token, ','); index = std::stoi(token);
            std::getline(ss, token, ','); name = token;  // satellite name
            std::getline(ss, token, ','); state.x = std::stod(token);
            std::getline(ss, token, ','); state.y = std::stod(token);
            std::getline(ss, token, ','); state.z = std::stod(token);
            std::getline(ss, token, ','); state.vx = std::stod(token);
            std::getline(ss, token, ','); state.vy = std::stod(token);
            std::getline(ss, token, ','); state.vz = std::stod(token);
            std::getline(ss, token, ','); state.orb.semi_major_axis = std::stod(token);
            std::getline(ss, token, ','); state.orb.eccentricity = std::stod(token);
            std::getline(ss, token, ','); state.orb.inclination = std::stod(token);
            std::getline(ss, token, ','); state.orb.RAAN = std::stod(token);
            std::getline(ss, token, ','); state.orb.augment_of_periapsis = std::stod(token);
            // true_anomaly is in the file but not needed for your function

            // Track time bounds
            min_time = std::min(min_time, state.time);
            max_time = std::max(max_time, state.time);

            // Store by satellite index
            sat_data[index].push_back(state);
        }

        file.close();

        std::cout << "Loaded data for " << sat_data.size() << " satellites" << std::endl;
        std::cout << "Time range: [" << min_time << ", " << max_time << "] seconds" << std::endl;
        std::cout << "Using tstep_size = " << tstep_size << " seconds (from Nbody.hpp)" << std::endl;
        std::cout << "Using t_final = " << t_final << " seconds (from Nbody.hpp)" << std::endl;
    }

    orbital_elements get_orbital_elements_at_time(int sat_idx, double time) const {
        const auto& trajectory = sat_data.at(sat_idx);

        // Binary search for time (assuming sorted by time)
        size_t left = 0;
        size_t right = trajectory.size() - 1;

        // Handle boundary cases
        if (time <= trajectory[0].time) {
            return trajectory[0].orb;
        }
        if (time >= trajectory[right].time) {
            return trajectory[right].orb;
        }

        // Find bracketing indices
        while (right - left > 1) {
            size_t mid = (left + right) / 2;
            if (trajectory[mid].time < time) {
                left = mid;
            } else {
                right = mid;
            }
        }

        // Linear interpolation between left and right
        double t0 = trajectory[left].time;
        double t1 = trajectory[right].time;
        double alpha = (time - t0) / (t1 - t0);

        orbital_elements orb;
        orb.semi_major_axis = (1 - alpha) * trajectory[left].orb.semi_major_axis +
                              alpha * trajectory[right].orb.semi_major_axis;
        orb.eccentricity = (1 - alpha) * trajectory[left].orb.eccentricity +
                          alpha * trajectory[right].orb.eccentricity;
        orb.inclination = (1 - alpha) * trajectory[left].orb.inclination +
                         alpha * trajectory[right].orb.inclination;

        // Angular elements need special care (wrap-around)
        orb.RAAN = interpolate_angle(trajectory[left].orb.RAAN,
                                     trajectory[right].orb.RAAN, alpha);
        orb.augment_of_periapsis = interpolate_angle(
            trajectory[left].orb.augment_of_periapsis,
            trajectory[right].orb.augment_of_periapsis, alpha);

        return orb;
    }

    // Helper to interpolate angles (handles 2π wrap-around)
    double interpolate_angle(double a0, double a1, double alpha) const {
        double diff = a1 - a0;

        // Wrap difference to [-π, π]
        while (diff > M_PI) diff -= 2 * M_PI;
        while (diff < -M_PI) diff += 2 * M_PI;

        double result = a0 + alpha * diff;

        // Wrap result to [0, 2π]
        while (result < 0) result += 2 * M_PI;
        while (result >= 2 * M_PI) result -= 2 * M_PI;

        return result;
    }
};
