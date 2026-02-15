cmake_minimum_required(VERSION 3.20)
project(Nbody)

set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc-15)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++-15)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Point CMake to Homebrew Armadillo
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/armadillo/15.0.3")
# Find Armadillo
find_package(Armadillo REQUIRED)

# Find Boost
find_package(Boost REQUIRED COMPONENTS math)

# Find CasADi using pkg-config
find_package(PkgConfig REQUIRED)
# Find OpenMP
find_package(OpenMP REQUIRED)


# Create a shared library for common code (excluding files with main functions)
set(LIBRARY_SOURCES
    "${PROJECT_SOURCE_DIR}/src/data_access_lib.cpp"
    "${PROJECT_SOURCE_DIR}/src/DP.cpp"
    "${PROJECT_SOURCE_DIR}/src/LambertSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/Local_search.cpp"
    "${PROJECT_SOURCE_DIR}/src/LowThrustAnalytical.cpp"
    "${PROJECT_SOURCE_DIR}/src/Nbody.cpp"
    "${PROJECT_SOURCE_DIR}/src/Simulated_Annealing.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tabu_search.cpp"
    "${PROJECT_SOURCE_DIR}/src/Trajectory_selection.cpp"
    "${PROJECT_SOURCE_DIR}/src/VNS.cpp"
)

add_library(NbodyLib STATIC ${LIBRARY_SOURCES})

# Include Armadillo headers for the library
target_include_directories(NbodyLib PRIVATE
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${Boost_INCLUDE_DIRS}
                            ${Casadi_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(NbodyLib PRIVATE ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${Casadi_LIBRARIES})

# Create separate executables
add_executable(ScalingTests "${PROJECT_SOURCE_DIR}/src/Sim.cpp")
add_executable(RunSimExecutable "${PROJECT_SOURCE_DIR}/src/RunSim.cpp")
#add_executable(SolverTrial "${PROJECT_SOURCE_DIR}/src/Solver_trial.cpp")

# Link executables to the library
target_include_directories(ScalingTests PRIVATE
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${Boost_INCLUDE_DIRS}
                            ${Casadi_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(ScalingTests PRIVATE NbodyLib ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${Casadi_LIBRARIES})

target_include_directories(RunSimExecutable PRIVATE
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${Boost_INCLUDE_DIRS}
                            ${Casadi_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
                            OpenMP::OpenMP_CXX
)
target_link_libraries(RunSimExecutable PRIVATE NbodyLib ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${Casadi_LIBRARIES})

#target_include_directories(SolverTrial PRIVATE
#                            ${ARMADILLO_INCLUDE_DIRS}
#                            ${Boost_INCLUDE_DIRS}
#                            ${Casadi_INCLUDE_DIRS}
#                            ${PROJECT_SOURCE_DIR}/include
#)
#target_link_libraries(SolverTrial PRIVATE NbodyLib ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES} ${Casadi_LIBRARIES})
