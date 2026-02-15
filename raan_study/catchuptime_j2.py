import matplotlib.pyplot as plt
import numpy as np 


def raanrate(a):
        
    return -3/2 * np.sqrt((3.986e14)/(6.378e6 + a)**3) * ((1.08263*10**(-3)) * np.cos(np.radians(56)))/(1-0.001**2)**2 * (6.378e6/ (6.378e6 + a))**2


def raan_phasing_time(a_target, a_phasing, delta_raan_initial):
    """
    Calculate time needed to close a RAAN gap using J2
    
    Parameters:
    - a_target: target orbit altitude (m)
    - a_phasing: phasing orbit altitude (m)
    - delta_raan_initial: initial RAAN difference (degrees)
    
    Returns: time in days
    """
    # RAAN rates
    rate_target = raanrate(a_target)  # rad/s
    rate_phasing = raanrate(a_phasing)  # rad/s
    
    # Differential rate
    diff_rate = rate_phasing - rate_target  # rad/s
    
    # Time to close gap
    time_seconds = np.radians(delta_raan_initial) / diff_rate
    
    return time_seconds / 86400


a_target = 700e3

a_phasing = np.linspace(800,1200)*1e3


rp = []

for ap in a_phasing:
    print(raan_phasing_time(a_target,ap,40))

    rp.append(raan_phasing_time(a_target,ap,40))
        

plt.plot(a_phasing/1e3,rp)

plt.xlabel("Depot Altitude (km)")
plt.ylabel("Phasing Time (days)")

plt.show()


