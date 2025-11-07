cmake_minimum_required(VERSION 3.20)
project(Nbody)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Point CMake to Homebrew Armadillo
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/armadillo/15.0.3")
# Find Armadillo
find_package(Armadillo REQUIRED)


file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCES})



# Include Armadillo headers and link library
target_include_directories(${PROJECT_NAME} PRIVATE 
                            ${ARMADILLO_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
)


target_link_libraries(${PROJECT_NAME} PRIVATE ${ARMADILLO_LIBRARIES})
