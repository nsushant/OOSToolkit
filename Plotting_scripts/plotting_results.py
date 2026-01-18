# visualize_satellites.py
# Usage: python3 plotting_scripts/visualize_satellites.py

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Folder of the current script (A/)
this_dir = os.path.dirname(os.path.abspath(__file__))

# Project root (parent of A/)
root = os.path.dirname(this_dir)

# Path to B/data.csv
data_path = os.path.join(root, "data", "WalkerDelta.csv")

print("Loading:", data_path)

df = pd.read_csv(data_path)

# Unique satellites 
sat_names = df['name'].unique()

# Setup 3D plot 
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Satellite Trajectories")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")

# Optional: plot Earth as a sphere for reference
import numpy as np
R_EARTH = 6378137.0
u = np.linspace(0, 2*np.pi, 50)
v = np.linspace(0, np.pi, 50)
x = R_EARTH * np.outer(np.cos(u), np.sin(v))
y = R_EARTH * np.outer(np.sin(u), np.sin(v))
z = R_EARTH * np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(x, y, z, color='b', alpha=0.2, zorder=0)

# Plot trajectories 
for name in sat_names:
    if (name == "service_1"):
        sat_data = df[df['name'] == name]
        ax.plot(sat_data['x'], sat_data['y'], sat_data['z'], '--')
    else:
        sat_data = df[df['name'] == name]
        ax.plot(sat_data['x'], sat_data['y'], sat_data['z'], label=name)

for name in sat_names:
    
    sat_data = df[df['name'] == name].iloc[0]
    ax.scatter3D(sat_data['x'], sat_data['y'], sat_data['z'])


# Optional: show legend (omit if too many satellites)
# ax.legend(loc='upper right', fontsize=8)
# Equal aspect ratio
ax.set_box_aspect([1,1,1])

plt.show()
