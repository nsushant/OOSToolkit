### On Orbit Servicing Toolkit

## Project Aims 

To create a simulation based toolkit that integrates ideas from OR and Physics to enable the study and planning of on-orbit servicing missions. The currently available components include ;- 

+ An N-body simulation that evolves the orbits of a walker delta constellation and a service depot forwards in time. 

+ Plotting scripts in python that allow you to visualize initial conditions, view orbits and animate the motions of constellation satellites 

+ A lambert solver that lets you optimize trajectories for service shuttles based on servicing demands. 

+ Optimization routines that let you construct both an optimal schedule and an optimal inventory management system.

## Installation

```
git clone https://github.com/nsushant/OOSToolkit.git
```

The package requires CMake, Armadillo (for vector calculations) and functions from the standard library. To get Armadillo you can follow the procedure at :

https://arma.sourceforge.net/download.html. 


## Usage 

The project uses CMake so if you would like to start with a fresh build, 

```
rm -rf build

mkdir build 

cd build

cmake ..

make

```

This will create an executable called Nbody which can be run from the build folder as follows. 

```
./Nbody 
```

Running this executable will trigger all functions within the Sim.cpp file and by default generate two output files that can be found in the /data directory. 

1) A file with position, velocity and time information about satellites simulated in Sim.cpp
2) A file containing all possible trajectories for a given schedule 


Plotting scripts written in python are also included which will allow you to visualize the fuel cost gradients, satellite trajectories and more. 


<div style="display: flex; justify-content: space-around;">
  <img src="./Figures/Figure_1.png" alt="trajectories image" width="45%">
  <img src="./Figures/DeltaV.jpeg" alt="DeltaV image" width="49%">
</div>

<p style="text-align: center; margin-top: 10px;">
    <strong>Figure 1:</strong> Example images created using the plotting script.The left panel shows a visualization of satellite trajectories while the right panel shows Delta V contours showing the fuel consumption for all possible trajectories between two satellites for a range of arrival and departure times. 
</p>
