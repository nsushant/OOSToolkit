# NbodyWalkerDelta C++ Library Documentation

## Overview

NbodyWalkerDelta is a comprehensive C++ library for satellite constellation simulation, trajectory optimization, and orbital mechanics calculations. The library provides tools for simulating Delta Walker satellite constellations, solving Lambert's problem for orbital transfers, and optimizing service schedules using various search algorithms.

## Dependencies

- **C++17** standard
- **Armadillo** linear algebra library (v15.0.3)
- CMake 3.20+

## Core Components

### 1. Nbody Simulation (`Nbody.hpp`)

Main module for satellite constellation dynamics and orbital mechanics.

#### Constants
- `MU_EARTH`: Earth's gravitational parameter (3.986004418e14 m³/s²)
- `R_EARTH`: Earth radius (6378137.0 m)
- `J2`: J2 perturbation coefficient (1.08262668e-3)
- `G_CONST`: Gravitational constant (6.67430e-11 m³/kg/s²)

#### Data Structures

**`orbital_elements`**
```cpp
struct orbital_elements {
    double inclination;          // Orbital inclination (rad)
    double RAAN;                 // Right Ascension of Ascending Node (rad)
    double augment_of_periapsis; // Argument of periapsis (rad)
    double true_anomaly;         // True anomaly (rad)
    double semi_major_axis;      // Semi-major axis (m)
    double eccentricity;          // Eccentricity
};
```

**`satellite_object`**
```cpp
struct satellite_object {
    std::string satname;  // Satellite identifier
    double mass;          // Mass (kg)
    arma::vec3 r;         // Position vector in ECI frame (m)
    arma::vec3 v;         // Velocity vector in ECI frame (m/s)
};
```

**`force_model`**
Immutable struct for configuring force model options. Once initialized, the boolean values cannot be changed.

```cpp
struct force_model {
    const bool includeJ2;       // Include J2 perturbations (immutable)
    const bool includeMutual;   // Include mutual gravitational forces (immutable)
    
    // Constructors
    force_model();                              // Default: both false
    force_model(bool j2, bool mutual);         // Custom values
};
```

#### Key Functions

- `deg_to_rads(double deg)`: Convert degrees to radians
- `get_inclination(arma::vec r, arma::vec v)`: Calculate orbital inclination
- `orb_elems_from_rv(arma::vec r, arma::vec v)`: Convert position/velocity to orbital elements
- `acceleration_due_to_central_body(const arma::vec3& r)`: Calculate central body acceleration
- `acceleration_due_to_J2(const arma::vec3& r)`: Calculate J2 perturbation acceleration
- `compute_acceleration(...)`: Compute total acceleration with force model options
- `runge_kutta_step(...)`: 4th-order Runge-Kutta integration step
- `build_walker_constellation(...)`: Create Walker delta constellation
- `run_simulation(...)`: Run complete orbital simulation

### 2. Lambert Solver (`LambertSolver.hpp`)

Implementation of Lambert's problem solver for orbital transfers, based on ESA's pykep library with corrections from D.Izzo's publication.

#### Key Functions

- `lambert_solver(...)`: Main solver returning velocity vectors for orbital transfer
  - **Inputs**: r1, r2 (position vectors), time_of_flight, mu, retrograde flag, max_revolutions
  - **Outputs**: Vector of velocity solutions
- `get_v1(...)`, `get_v2(...)`: Extract velocity vectors from solutions
- `taking_derivatives(...)`: Compute derivatives for Householder iteration
- `tof_lagrange(...)`: Calculate time of flight using Lagrange coefficients
- `hypergeometricF(...)`: Hypergeometric function calculation
- `iterate_householder(...)`: Householder iteration for convergence

### 3. Low-Thrust Analytical (`LowThrustAnalytical.hpp`)

Analytical methods for low-thrust trajectory optimization.

#### Key Functions

- `calculate_plane_diff_angle(...)`: Calculate angle difference between orbital planes
- `calculate_edelbaum_deltaV(...)`: Calculate Edelbaum delta-V for low-thrust transfers
- `low_thrust_di_carlo_independent()`: Independent Monte Carlo analysis
- `low_thrust_di_carlo_twophase()`: Two-phase Monte Carlo analysis

### 4. Trajectory Selection (`Trajectory_selection.hpp`)

Optimal trajectory selection algorithms for satellite servicing missions.

#### Key Functions

- `find_optimal_trajectory(...)`: Find optimal Lambert parameters for minimum delta-V
- `compact_optimal_calc(...)`: Compact optimal calculation between two satellites
- `find_optimal_trajectory_no_iter(...)`: Non-iterative optimal trajectory finding
- `run_exhaustive_search(...)`: Exhaustive search for optimal trajectories

### 5. Local Search (`Local_search.hpp`)

Local search algorithms for schedule optimization.

#### Data Structures

**`task_block`**
```cpp
struct task_block {
    double arrival_time;
    double departure_time;
    std::string satname;
    double deltaV_arrival;
    double service_duration;
    double arrival_constraint;
};
```

**`schedule_struct`**
```cpp
struct schedule_struct {
    std::vector<task_block> blocks;
};
```

#### Key Functions

- `move_dt(...)`, `move_dt2(...)`: Time adjustment moves
- `move_add_arrival(...)`, `move_add_departure(...)`: Arrival/departure adjustment moves
- `swap_slots(...)`: Slot swapping operations
- `create_schedule(...)`, `create_schedule_lambert_only(...)`: Schedule creation
- `local_search_opt_schedule_lambert_only(...)`: Local search optimization
- `run_local_search(...)`: Execute local search algorithm
- `init_dep_arrival_times_random(...)`: Random initialization
- `init_dep_arrival_times_strict_timespan(...)`: Strict time span initialization

### 6. Variable Neighborhood Search (`VNS.hpp`)

Variable Neighborhood Search algorithm for global optimization.

#### Key Functions

- `vn_search(...)`: Variable neighborhood search implementation
- `run_vn_search(...)`: Execute VNS with flexible arrival times
- `run_vn_search_fixed_tarrive(...)`: Execute VNS with fixed arrival times

### 7. Tabu Search (`Tabu_search.hpp`)

Tabu search algorithm for combinatorial optimization.

#### Key Functions

- `Tabu_search_opt_schedule_lambert_only(...)`: Tabu search optimization
- `run_tabu_search(...)`: Execute tabu search algorithm

### 8. Simulated Annealing (`Simulated_Annealing.hpp`)

Simulated annealing optimization algorithm.

#### Key Functions

- `prob_calculation(...)`: Calculate acceptance probability
- `simulated_annealing_lambert(...)`: Simulated annealing with Lambert transfers
- `run_simulated_annealing(...)`: Execute simulated annealing algorithm

### 9. Exact Methods (`Exact_methods.hpp`)

Exact optimization methods using dynamic programming.

#### Key Functions

- `is_feasible_sol(...)`: Check solution feasibility
- `finding_individual_minimas_dynamic_programming(...)`: Dynamic programming for minima
- `calculate_deltav_upto_thisblock(...)`: Calculate cumulative delta-V
- `dynamic_program_fixed_tasksize_Tfixed(...)`: Fixed task size dynamic programming

### 10. Data Access Library (`data_access_lib.hpp`)

CSV data handling with pandas-like interface.

#### Data Structure

**`DataFrame`**
```cpp
struct DataFrame {
    std::unordered_map<std::string, std::vector<std::string>> data;
    std::vector<std::string> headers;
};
```

#### Key Functions

- `DataFrame(const std::string &filename)`: Constructor from CSV file
- `operator[](const std::string &colname)`: Column access (pandas-style)
- `getNumeric(const std::string &colname)`: Convert column to numeric vector
- `showColumns()`: Display available columns
- `find_idxs_of_match(...)`: Find indices of matching values

## Usage Examples

### Building the Library

```bash
mkdir build && cd build
cmake ..
make
```

### Basic Simulation

```cpp
#include "Nbody.hpp"

// Create Walker constellation
std::vector<satellite_object> constellation = build_walker_constellation(
    3,      // num_planes
    24,     // total_satnum
    1,      // phase
    500000, // altitude_m
    deg_to_rads(53), // inclination_rad
    100     // sat_mass
);

// Run simulation
run_simulation("output.csv", "walker_delta", 86400, 10, 500000, 3, 24, 1, deg_to_rads(53));
```

### Lambert Transfer

```cpp
#include "LambertSolver.hpp"

arma::vec3 r1 = {7000e3, 0, 0};  // Initial position
arma::vec3 r2 = {0, 7000e3, 0};  // Final position
double tof = 3600;               // Time of flight (s)

std::vector<arma::vec> solutions = lambert_solver(r1, r2, tof, MU_EARTH, 1, 0);
arma::vec3 v1 = get_v1(solutions, 0);  // Initial velocity
arma::vec3 v2 = get_v2(solutions, 0);  // Final velocity
```

### Schedule Optimization

```cpp
#include "Local_search.hpp"
#include "data_access_lib.hpp"

// Load simulation data
DataFrame simfile("simulation_data.csv");

// Create initial schedule
schedule_struct schedule = create_schedule_lambert_only(
    deltaV_of_schedule_init, arrival_times, departure_times, satnames, simfile
);

// Optimize with local search
schedule_struct optimized = local_search_opt_schedule_lambert_only(
    init_deltaV, schedule, dt_move, simfile, service_time, move_methods
);
```

## File Structure

```
NbodyWalkerDelta/
├── include/
│   ├── Nbody.hpp
│   ├── LambertSolver.hpp
│   ├── LowThrustAnalytical.hpp
│   ├── Trajectory_selection.hpp
│   ├── Local_search.hpp
│   ├── VNS.hpp
│   ├── Tabu_search.hpp
│   ├── Simulated_Annealing.hpp
│   ├── Exact_methods.hpp
│   └── data_access_lib.hpp
├── src/
│   ├── Nbody.cpp
│   ├── LambertSolver.cpp
│   ├── LowThrustAnalytical.cpp
│   ├── Trajectory_selection.cpp
│   ├── Local_search.cpp
│   ├── VNS.cpp
│   ├── Tabu_search.cpp
│   ├── Simulated_Annealing.cpp
│   ├── Exact_methods.cpp
│   └── data_access_lib.cpp
├── CMakeLists.txt
└── docs.md (this file)
```

## Units and Conventions

- **Units**: meters (m), seconds (s), kilograms (kg)
- **Angles**: Radians (use `deg_to_rads()` for conversion)
- **Coordinate System**: ECI (Earth-Centered Inertial)
- **Time Integration**: 4th-order Runge-Kutta

## Author

S. Nigudkar (2025)

## References

1. D. Izzo, "Revisiting Lambert's Problem", ESA ACT-MAD (2014)
2. pykep library: https://github.com/esa/pykep
3. ESA Advanced Concepts Team