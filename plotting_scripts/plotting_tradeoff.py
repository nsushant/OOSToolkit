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

mpl.style.use('dark_background')

# Under Construction


def calculate_nondim_T(r1 , r2, mu, tof):

    c = r2 - r1

    r1_norm = np.sqrt(r1[0]**2 + r1[1]**2 + r1[2]**2)
    r2_norm = np.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)

    c = np.sqrt(c[0]**2 + c[1]**2 + c[2]**2)

    s = 1/2 * (r1_norm + r2_norm + c)
    l = np.sqrt(1 - c/s)

    T = np.sqrt(2*mu / s**3) * tof
    T_1 = 2/3 * (1 - l**3)
    T_00 = np.acos(l) + l*np.sqrt(1-l)

    return T, T_1, T_00, l



# Folder of the current script (A/)
this_dir = os.path.dirname(os.path.abspath(__file__))

# Project root (parent of A/)
root = os.path.dirname(this_dir)

# Path to B/data.csv
data_path = os.path.join(root, "data", "Trajectories.csv")

csvfile = pd.read_csv(data_path)

minima = csvfile[csvfile['deltaV'] == csvfile['deltaV'].min()]

tdep = csvfile['t_depot'].values
tarrive = csvfile['t_client'].values

tdep_grid = np.linspace(min(tdep), max(tdep), 100)
tarrive_grid = np.linspace(min(tarrive), max(tarrive), 100)
T_DEP, T_ARR = np.meshgrid(tdep_grid, tarrive_grid)


deltaV_grid = griddata(
    (tdep, tarrive),   # input points
    csvfile['deltaV'].values,            # values
    (T_DEP, T_ARR),    # output grid
    method='cubic'     # or 'linear'
)

mu = 3.986004418e14

Ts = []
T1s = []
T00s = []
ls = []
for i in range(len(csvfile['x_ser'].values)): 
    row = csvfile.iloc[i]
    r1 = np.array([row['x_ser'],row['y_ser'],row['z_ser']])
    r2 = np.array([row['x_cl'],row['y_cl'],row['z_cl']])
    tof = row['t_client'] - row['t_depot']

    T, T_1, T_00, l = calculate_nondim_T(r1 , r2, mu,tof)
    Ts.append(T)
    T1s.append(T_1)
    T00s.append(T_00)
    ls.append(l)

'''

plt.scatter(tdep,tarrive,c=ls)
plt.colorbar(label="$\lambda$")
plt.xlabel('Departure Time [seconds]')
plt.ylabel('Arrival Time [seconds]')
plt.scatter(1200,5000,c='yellow',label="optimal transfer")
plt.legend()

'''


T_grid = griddata(
    (tdep, tarrive),   # input points
    np.asarray(Ts), # values
    (T_DEP, T_ARR),    # output grid
    method='cubic'     # or 'linear'
)

T1_grid = griddata(
    (tdep, tarrive),   # input points
    np.asarray(T1s), # values
    (T_DEP, T_ARR),    # output grid
    method='cubic'     # or 'linear'
)

T0_grid = griddata(
    (tdep, tarrive),   # input points
    np.asarray(T00s), # values
    (T_DEP, T_ARR),    # output grid
    method='cubic'     # or 'linear'
)

l_grid = griddata(
    (tdep, tarrive),   # input points
    np.asarray(ls), # values
    (T_DEP, T_ARR),    # output grid
    method='cubic'     # or 'linear'
)
'''

rdiff = (np.sqrt((csvfile['x_ser']-csvfile['x_cl'])**2+(csvfile['y_ser']-csvfile['y_cl'])**2+(csvfile['z_ser']-csvfile['z_ser'])**2))*10**(-6)

plt.plot(csvfile['t_client'],rdiff/max(rdiff),label="Travel Distance")

plt.plot(csvfile['t_client'],csvfile['deltaV']/csvfile['deltaV'].max(),label="$\Delta V$")

plt.xlabel("T_client")

plt.title("Normalized Params when T_depot = 0")

plt.legend()
'''







plt.figure(figsize=(10,6))

vmin, vmax = np.nanmin(deltaV_grid), np.nanmax(deltaV_grid)
levels = np.logspace(np.log10(vmin), np.log10(vmax), 40)

#levels = np.linspace(np.nanmin(deltaV_grid), np.nanmax(deltaV_grid), 40)
contour = plt.contourf(T_DEP, T_ARR, deltaV_grid, levels=levels, cmap='plasma',norm=colors.LogNorm(vmin=vmin, vmax=vmax))

#cs = plt.contour(T_DEP, T_ARR, deltaV_grid, colors='black', linewidths=1.0,norm=colors.LogNorm(vmin=vmin, vmax=vmax))
cs = plt.contour(T_DEP, T_ARR, l_grid,colors='black')

plt.clabel(cs, inline=True, fontsize=8, fmt='%.1f')
plt.xlabel('Departure Time [seconds]')
plt.ylabel('Arrival Time [seconds]')
plt.title('Transfer energies over one orbit (6000s)')
plt.scatter(1200,5000,c='yellow',label="optimal transfer")
plt.scatter(minima['t_depot'].values,minima['t_client'].values,marker="+",c='r',label='Actual Minima')
cbar = plt.colorbar(contour, label='ΔV [km/s]',norm = colors.LogNorm(vmin=vmin, vmax=vmax))
cbar.locator = LogLocator(base=10)    # sets ticks at 10^n
cbar.update_ticks()
plt.legend()
plt.grid(True, linestyle=':')

plt.tight_layout()

#plt.scatter(csvfile['tof'],csvfile['deltaV'])



plt.show()