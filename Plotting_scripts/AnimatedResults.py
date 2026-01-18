import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting
import os 
import argparse
import sys

def load_data(csv_path):
    """Load data from CSV file"""
    print("Loading:", csv_path)
    try:
        df = pd.read_csv(csv_path)
        return df
    except Exception as e:
        print(f"Error loading {csv_path}: {e}")
        sys.exit(1)

# Get data path from command line or use default
if len(sys.argv) > 1:
    csv_path = sys.argv[1]
else:
    # Folder of the current script 
    this_dir = os.path.dirname(os.path.abspath(__file__))
    # Project root 
    root = os.path.dirname(this_dir)
    # Path to data 
    csv_path = os.path.join(root, "data", "WalkerDelta.csv")

df = load_data(csv_path)

# Get satellite info
satellites = df['name'].unique()
times = sorted(df['time_s'].unique())
n_steps = len(times)

# plot figure
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')

# Plot Earth as reference sphere
R_earth = 6378137.0
u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
x = R_earth*np.cos(u)*np.sin(v)
y = R_earth*np.sin(u)*np.sin(v)
z = R_earth*np.cos(v)
ax.plot_surface(x, y, z, color='lightblue', alpha=0.3, linewidth=0)

# Set equal aspect ratio
max_radius = df[['x','y','z']].abs().values.max()
for axis in 'xyz':
    getattr(ax, f'set_{axis}lim')([-max_radius, max_radius])

# Create a scatter object for satellites
scatter = ax.scatter([], [], [], s=20, c='red')


for name in satellites:
    if (name == "service_1"):
        sat_data = df[df['name'] == name]
        ax.plot(sat_data['x'], sat_data['y'], sat_data['z'], '--')
    else:
        sat_data = df[df['name'] == name]
        ax.plot(sat_data['x'], sat_data['y'], sat_data['z'], label=name)

# === Animation Update Function ===
def update(frame):
    t = times[frame]
    current = df[df['time_s'] == t]
    scatter._offsets3d = (current['x'], current['y'], current['z'])
    ax.set_title(f"Time = {t:.0f} s")
    return scatter,

# === Animate ===
ani = FuncAnimation(fig, update, frames=n_steps, interval=100, blit=False)

plt.show()
