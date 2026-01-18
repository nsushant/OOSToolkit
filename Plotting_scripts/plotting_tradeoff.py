import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.stats import gaussian_kde
import os
import numpy as np 
#import seaborn as sns
from scipy.interpolate import griddata
import matplotlib.colors as colors
from matplotlib.ticker import LogLocator
import matplotlib as mpl

mpl.style.use('default')

# Under Construction


# Folder of the current script (A/)
this_dir = os.path.dirname(os.path.abspath(__file__))

# Project root (parent of A/)
root = os.path.dirname(this_dir)

# Path to B/data.csv
data_path = os.path.join(root, "data", "DeltaV_vs_movesize.csv")

csvfile = pd.read_csv(data_path)

print(csvfile.columns)
plt.plot(csvfile["move_size"], csvfile['deltaV_improvement']/csvfile['deltaV_improvement'].max(),label= "Delta V",lw=2)         
plt.plot(csvfile["move_size"], csvfile['time_improvement']/csvfile['time_improvement'].max(),label= "Timespan",lw=2)

plt.title("Move : Earlier Arrival Times, Later Departures")
plt.xlabel("Move size (s)")
plt.ylabel("Improvements")
plt.legend(frameon=False)
plt.show()