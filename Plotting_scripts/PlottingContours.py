import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
import os 
import sys 
import matplotlib as mpl 


mpl.style.use('dark_background')


data_path = os.path.join("/Users/sushantnigudkar/Downloads/NbodyWalkerDelta", "data", "Trajectories.csv")
csvfile = pd.read_csv(data_path)
csvfile = csvfile[csvfile['t_client']<6000]


departure_times = csvfile['t_depot']
arrival_times = csvfile['t_client']
fuel_costs = csvfile['deltaV']

minima = np.argmin(fuel_costs.values)


depmin = departure_times.values[minima]
arrmin = arrival_times.values[minima]


# Create a regular grid of points
dep_grid, arr_grid = np.meshgrid(
    np.linspace(departure_times.min(), departure_times.max(), 100),
    np.linspace(arrival_times.min(), arrival_times.max(), 100)
)

# Interpolate your 1D data onto the grid
fuel_grid = griddata(
    (departure_times, arrival_times),
    fuel_costs,
    (dep_grid, arr_grid),
    method='cubic'
)

# Define log-spaced levels
levels = np.logspace(
    np.log10(fuel_costs.min()),
    np.log10(fuel_costs.max()),
    num=20
)

# Plot contours with log scale
plt.figure(figsize=(8, 6))
contour = plt.contourf(
    dep_grid, arr_grid, fuel_grid,
    levels=levels, cmap='inferno', norm=LogNorm()
)

#  Add colorbar and labels
cbar = plt.colorbar(contour,norm=LogNorm,ticks=[1e3,1e4,1e5,1e6])
cbar.set_label('Fuel Cost (log scale)')


plt.scatter(0.0,4800.0)

plt.xlabel('Departure Time From Depot  (s)')
plt.ylabel('Arrival Time At Client (s)')
plt.title('$\Delta V$ Contours ')
plt.tight_layout()
plt.savefig("Deltav_contours.pdf")

