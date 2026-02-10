#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <stdexcept>

// Physical constants
constexpr double MU_EARTH = 398600.4418;        // km^3/s^2
constexpr double R_EARTH = 6378.137;            // km
constexpr double J2 = 0.00108263;               // J2 coefficient
constexpr double PI = 3.14159265358979323846;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;

// Orbital elements structure
struct OrbitalElements {
    double a;       // Semi-major axis (km)
    double i;       // Inclination (rad)
    double Omega;   // RAAN (rad)

    OrbitalElements(double _a, double _i, double _Omega)
        : a(_a), i(_i), Omega(_Omega) {}
};

// Solution structure
struct DiCarloSolution {
    double psi_bar;      // Thrust arc semi-amplitude (rad)
    double beta_bar;     // Elevation angle (rad)
    double ToF1;         // Phase 1 time of flight (days)
    double ToF2;         // Phase 2 time of flight (days)
    double delta_v;      // Total delta-V (km/s)
    bool converged;      // Did the solver converge?

    void print() const {
        std::cout << "=== Di Carlo Strategy 1 Solution ===" << std::endl;
        std::cout << "Psi (thrust arc):  " << psi_bar * RAD_TO_DEG << " deg" << std::endl;
        std::cout << "Beta (elevation):  " << beta_bar * RAD_TO_DEG << " deg" << std::endl;
        std::cout << "Phase 1 (coast):   " << ToF1 << " days" << std::endl;
        std::cout << "Phase 2 (thrust):  " << ToF2 << " days" << std::endl;
        std::cout << "Total Delta-V:     " << delta_v << " km/s" << std::endl;
        std::cout << "Converged:         " << (converged ? "Yes" : "No") << std::endl;
    }
};

// Simple Newton-Raphson solver for 1D root finding
class NewtonRaphson {
public:

    //function exists for the whole duration of the program
    static double solve(std::function<double(double)> f,
                       std::function<double(double)> df,
                       double x0,
                       double tol = 1e-9,
                       int max_iter = 100) {

        // the aim is to find minima/root f(x) = 0;
        double x = x0; // set initial guess
        for (int i = 0; i < max_iter; ++i) {
            double fx = f(x); //eval function at init guess

            if (std::abs(fx) < tol) {

                return x;

            }

            double dfx = df(x); // take derivative

            if (std::abs(dfx) < 1e-12) {
                throw std::runtime_error("Newton-Raphson: derivative too small");
            }
            x = x - fx / dfx; // this is the correction that was the output of our current iteration
        }
        throw std::runtime_error("Newton-Raphson: failed to converge");
    }
};



// Di Carlo Strategy 1 implementation
class DiCarloStrategy1 {

private:
    OrbitalElements initial_;
    OrbitalElements final_;
    double ToF_total_;  // days
    double epsilon_;    // m/s^2

    // Helper functions
    double compute_raan_rate(double a, double i) const {
        // Secular RAAN rate due to J2 (rad/s)
        double n = std::sqrt(MU_EARTH / (a * a * a));
        return -1.5 * n * J2 * std::pow(R_EARTH / a, 2) * std::cos(i);
    }

    double compute_k1(double psi_bar, double beta_bar) const {
        double a0 = initial_.a;
        double af = final_.a;
        double i0 = initial_.i;
        double if_ = final_.i;

        double log_ratio = std::log(af / a0);
        double delta_i = if_ - i0;

        if (std::abs(delta_i) < 1e-10) {
            // Special case: no inclination change (Eq. 90)
            return 0.0;
        }

        double exp_term = std::exp(4.0 * log_ratio * i0 / delta_i);

        return -(3.0 * PI * MU_EARTH * J2 * R_EARTH * R_EARTH) /
               (4.0 * std::pow(a0, 4) * epsilon_ * std::sin(beta_bar) *
                std::sin(psi_bar)) * exp_term;
    }

    double compute_k2() const {
        double a0 = initial_.a;
        double af = final_.a;
        double i0 = initial_.i;
        double if_ = final_.i;

        double log_ratio = std::log(af / a0);
        double delta_i = if_ - i0;

        if (std::abs(delta_i) < 1e-10) {
            return 0.0;
        }

        return -4.0 * log_ratio / delta_i;
    }

    double compute_k3(double k2) const {
        double i0 = initial_.i;
        double if_ = final_.i;

        double term0 = std::exp(k2 * i0) * (k2 * std::cos(i0) + std::sin(i0));
        double term_f = std::exp(k2 * if_) * (k2 * std::cos(if_) + std::sin(if_));

        return term_f - term0;
    }

    double compute_Omega_2f(double psi_bar, double beta_bar, double ToF1) const {
        // RAAN at end of phase 2 (Eq. 88)
        double Omega_0 = initial_.Omega;
        double i0 = initial_.i;
        double a0 = initial_.a;

        // Phase 1 drift
        double Omega_drift_phase1 = -1.5 * std::sqrt(MU_EARTH) * J2 *
                                    R_EARTH * R_EARTH * std::cos(i0) *
                                    std::pow(a0, -3.5) * ToF1 * 86400.0;  // Convert days to seconds

        double Omega_1f = Omega_0 + Omega_drift_phase1;

        // Phase 2 contribution from J2
        double k1 = compute_k1(psi_bar, beta_bar);
        double k2 = compute_k2();
        double k3 = compute_k3(k2);

        if (std::abs(final_.i - initial_.i) < 1e-10) {
            // Special case: no inclination change (Eq. 90)
            double a0 = initial_.a;
            double af = final_.a;
            double i = initial_.i;

            double Omega_phase2 = (3.0 * PI * MU_EARTH * J2 * R_EARTH * R_EARTH * std::cos(i)) /
                                 (32.0 * epsilon_ * psi_bar) *
                                 (1.0 / std::pow(af, 4) - 1.0 / std::pow(a0, 4));

            return Omega_1f + Omega_phase2;
        }

        double Omega_phase2 = k1 / (1.0 + k2 * k2) * k3;

        return Omega_1f + Omega_phase2;
    }

    double compute_ToF2(double psi_bar, double beta_bar) const {
        // Time of flight for phase 2 (Eq. 84), returns seconds
        double a0 = initial_.a;
        double af = final_.a;
        double i0 = initial_.i;
        double if_ = final_.i;

        double log_ratio = std::log(af / a0);
        double delta_i = if_ - i0;

        double sqrt_term = 1.0;
        if (std::abs(delta_i) > 1e-10 && std::abs(log_ratio) > 1e-10) {
            sqrt_term = std::sqrt(1.0 + (4.0 * psi_bar * psi_bar * delta_i * delta_i) /
                                  (std::sin(psi_bar) * std::sin(psi_bar) * log_ratio * log_ratio));
        }

        double ToF2_sec = (PI / (2.0 * psi_bar)) *
                         (std::sqrt(MU_EARTH / a0) - std::sqrt(MU_EARTH / af)) *
                         sqrt_term;

        return ToF2_sec / 86400.0;  // Convert to days
    }

    double compute_beta_from_psi(double psi_bar) const {
        // Eq. 82: Relationship between beta and psi
        double a0 = initial_.a;
        double af = final_.a;
        double i0 = initial_.i;
        double if_ = final_.i;

        double log_ratio = std::log(af / a0);
        double delta_i = if_ - i0;

        if (std::abs(log_ratio) < 1e-10) {
            // No semi-major axis change
            return PI / 2.0;  // Pure out-of-plane
        }

        if (std::abs(delta_i) < 1e-10) {
            // No inclination change
            return 0.0;  // Pure in-plane
        }

        double tan_beta = (2.0 * psi_bar * delta_i) / (std::sin(psi_bar) * log_ratio);
        return std::atan(tan_beta);
    }

public:
    DiCarloStrategy1(const OrbitalElements& initial,
                     const OrbitalElements& final,
                     double ToF_days,
                     double epsilon_ms2)
        : initial_(initial), final_(final),
          ToF_total_(ToF_days), epsilon_(epsilon_ms2 * 1e-3) {}  // Convert to km/s^2

    DiCarloSolution solve() {
        DiCarloSolution solution;
        solution.converged = false;

        try {
            // Define the equation to solve: ToF1 + ToF2 = ToF_total with Omega_f constraint
            // We'll solve for psi_bar that satisfies both constraints

            auto objective = [this](double psi) -> double {
                // Given psi, compute beta
                double beta = compute_beta_from_psi(psi);

                // Compute ToF2
                double ToF2 = compute_ToF2(psi, beta);

                // Compute ToF1 from total time constraint
                double ToF1 = ToF_total_ - ToF2;

                if (ToF1 < 0) {
                    return 1e10;  // Infeasible
                }

                // Compute final RAAN
                double Omega_computed = compute_Omega_2f(psi, beta, ToF1);

                // Return error in RAAN
                return Omega_computed - final_.Omega;
            };

            // Numerical derivative
            auto derivative = [&objective](double psi) -> double {
                double h = 1e-6;
                return (objective(psi + h) - objective(psi - h)) / (2.0 * h);
            };

            // Initial guess for psi_bar
            double psi_initial = PI / 4.0;  // 45 degrees

            // Solve for psi_bar
            solution.psi_bar = NewtonRaphson::solve(objective, derivative, psi_initial);

            // Compute corresponding parameters
            solution.beta_bar = compute_beta_from_psi(solution.psi_bar);
            solution.ToF2 = compute_ToF2(solution.psi_bar, solution.beta_bar);
            solution.ToF1 = ToF_total_ - solution.ToF2;

            // Compute delta-V (Eq. 93)
            double a0 = initial_.a;
            double af = final_.a;
            double i0 = initial_.i;
            double if_ = final_.i;

            double log_ratio = std::log(af / a0);
            double delta_i = if_ - i0;

            double sqrt_term = 1.0;
            if (std::abs(delta_i) > 1e-10 && std::abs(log_ratio) > 1e-10) {
                sqrt_term = std::sqrt(1.0 +
                    (4.0 * solution.psi_bar * solution.psi_bar * delta_i * delta_i) /
                    (std::sin(solution.psi_bar) * std::sin(solution.psi_bar) *
                     log_ratio * log_ratio));
            }

            solution.delta_v = (std::sqrt(MU_EARTH / a0) - std::sqrt(MU_EARTH / af)) *
                              sqrt_term;

            solution.converged = true;

        } catch (const std::exception& e) {
            std::cerr << "Solver failed: " << e.what() << std::endl;
            solution.converged = false;
        }

        return solution;
    }

    // Compute trajectory evolution
    std::vector<OrbitalElements> compute_trajectory(const DiCarloSolution& sol,
                                                    int num_points = 100) const {
        std::vector<OrbitalElements> trajectory;

        // Phase 1: Coast (J2 drift only)
        int points_phase1 = static_cast<int>(num_points * sol.ToF1 / ToF_total_);
        for (int i = 0; i <= points_phase1; ++i) {
            double t = sol.ToF1 * i / points_phase1;  // days
            double t_sec = t * 86400.0;  // seconds

            // RAAN drift due to J2 (Eq. 77)
            double Omega = initial_.Omega - 1.5 * std::sqrt(MU_EARTH) * J2 *
                          R_EARTH * R_EARTH * std::cos(initial_.i) *
                          std::pow(initial_.a, -3.5) * t_sec;

            trajectory.emplace_back(initial_.a, initial_.i, Omega);
        }

        // Phase 2: Thrust (Eqs. 83, 87)
        int points_phase2 = num_points - points_phase1;
        for (int i = 1; i <= points_phase2; ++i) {
            double t = sol.ToF2 * i / points_phase2;  // days from start of phase 2
            double t_sec = t * 86400.0;  // seconds

            // Semi-major axis evolution (Eq. 83)
            double term1 = 1.0 / std::sqrt(initial_.a);
            double term2 = (2.0 * std::cos(sol.beta_bar) * sol.psi_bar * epsilon_ * t_sec) / PI;
            double denominator = 1.0 - 2.0 * std::sqrt(initial_.a / MU_EARTH) * term2;

            double a = MU_EARTH * std::pow(term1 + term2, 2) * std::pow(denominator, -2);

            // Inclination evolution (Eq. 83)
            double log_a_ratio = std::log(a / initial_.a);
            double i = initial_.i + (std::tan(sol.beta_bar) * std::sin(sol.psi_bar) *
                                    log_a_ratio) / (2.0 * sol.psi_bar);

            // RAAN evolution (integrate Eq. 85)
            // Simplified: use average a and i for this segment
            double Omega_prev = trajectory.back().Omega;
            double dt_sec = (sol.ToF2 / points_phase2) * 86400.0;
            double dOmega = -1.5 * std::sqrt(MU_EARTH) * J2 * R_EARTH * R_EARTH *
                           std::cos(i) * std::pow(a, -3.5) * dt_sec;

            trajectory.emplace_back(a, i, Omega_prev + dOmega);
        }

        return trajectory;
    }
};

// Example usage
int main() {
    try {
        // Example: LEO debris removal with large RAAN change
        // Initial orbit: 800 km altitude, 98 deg inclination, 0 deg RAAN
        // Final orbit: 900 km altitude, 99 deg inclination, 150 deg RAAN

        double a0 = 6378.137 + 800.0;  // km
        double af = 6378.137 + 900.0;  // km
        double i0 = 98.0 * DEG_TO_RAD;
        double if_ = 99.0 * DEG_TO_RAD;
        double Omega0 = 0.0 * DEG_TO_RAD;
        double Omegaf = 150.0 * DEG_TO_RAD;

        OrbitalElements initial(a0, i0, Omega0);
        OrbitalElements final(af, if_, Omegaf);

        double ToF_days = 600.0;
        double epsilon_ms2 = 1.5e-4;  // m/s^2

        std::cout << "=== Transfer Parameters ===" << std::endl;
        std::cout << "Initial: a=" << a0 << " km, i=" << i0*RAD_TO_DEG
                  << " deg, Omega=" << Omega0*RAD_TO_DEG << " deg" << std::endl;
        std::cout << "Final:   a=" << af << " km, i=" << if_*RAD_TO_DEG
                  << " deg, Omega=" << Omegaf*RAD_TO_DEG << " deg" << std::endl;
        std::cout << "Time of Flight: " << ToF_days << " days" << std::endl;
        std::cout << "Thrust: " << epsilon_ms2 << " m/s^2" << std::endl;
        std::cout << std::endl;

        DiCarloStrategy1 solver(initial, final, ToF_days, epsilon_ms2);
        DiCarloSolution solution = solver.solve();

        solution.print();

        if (solution.converged) {
            std::cout << "\n=== Generating Trajectory ===" << std::endl;
            auto trajectory = solver.compute_trajectory(solution, 20);

            std::cout << "Time(days)\ta(km)\t\ti(deg)\t\tOmega(deg)" << std::endl;
            for (size_t i = 0; i < trajectory.size(); ++i) {
                double t = ToF_days * i / (trajectory.size() - 1);
                std::cout << t << "\t\t"
                         << trajectory[i].a << "\t"
                         << trajectory[i].i * RAD_TO_DEG << "\t\t"
                         << trajectory[i].Omega * RAD_TO_DEG << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
