cmake_minimum_required(VERSION 3.20)
project(Nbody)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Point CMake to Homebrew Armadillo
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/armadillo/15.0.3")
# Find Armadillo
find_package(Armadillo REQUIRED)


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
                            ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(NbodyLib PRIVATE ${ARMADILLO_LIBRARIES})

# Create separate executables
add_executable(SimExecutable "${PROJECT_SOURCE_DIR}/src/Sim.cpp")
add_executable(RunSimExecutable "${PROJECT_SOURCE_DIR}/src/RunSim.cpp")

# Link executables to the library
target_include_directories(SimExecutable PRIVATE
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(SimExecutable PRIVATE NbodyLib ${ARMADILLO_LIBRARIES})

target_include_directories(RunSimExecutable PRIVATE
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(RunSimExecutable PRIVATE NbodyLib ${ARMADILLO_LIBRARIES})
