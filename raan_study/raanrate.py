import numpy as np 
import matplotlib.pyplot as plt 

def raanrate(a):
        
    return -3/2 * np.sqrt((3.986e14)/(6.378e6 + a)**3) * ((1.08263*10**(-3)) * np.cos(np.radians(56)))/(1-0.001**2)**2 * (6.378e6/ (6.378e6 + a))**2



t = 86400 

sched = np.arange(1,15)

ts = sched*t 

plt.plot(sched,abs(np.degrees(raanrate(600e3)*ts)))
plt.text(sched[-1], abs(np.degrees(raanrate(600e3)*ts))[-1], f'  {600e3/1e3:.0f} km')  #, va='center')

plt.plot(sched,abs(np.degrees(raanrate(800e3)*ts)))
plt.text(sched[-1], abs(np.degrees(raanrate(800e3)*ts))[-1], f'  {800e3/1e3:.0f} km')  #, va='center')

plt.plot(sched,abs(np.degrees(raanrate(1000e3)*ts)))
plt.text(sched[-1], abs(np.degrees(raanrate(1000e3)*ts))[-1], f'  {1000e3/1e3:.0f} km')  #, va='center')

plt.plot(sched,abs(np.degrees(raanrate(1200e3)*ts)))
plt.text(sched[-1], abs(np.degrees(raanrate(1200e3)*ts))[-1], f'  {1200e3/1e3:.0f} km')  #, va='center')

plt.xlim(0,18)


plt.title("$\Omega_{J2}$ Vs. Altitude")
plt.ylabel("$\Delta$ RAAN (deg)")
plt.xlabel("Days")
plt.show()
 
