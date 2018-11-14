import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%%
def func_10(E, V_0):
    theta_p = np.deg2rad(10)
    return np.rad2deg(np.arcsin(np.sqrt((E + V_0)/E) * np.sin(theta_p)) - theta_p)

def func_30(E, V_0):
    theta_p = np.deg2rad(30)
    return np.rad2deg(np.arcsin(np.sqrt((E + V_0)/E) * np.sin(theta_p)) - theta_p)
    
#%%
theta_10 = np.loadtxt('theta_10.csv')
theta_30 = np.loadtxt('theta_30.csv')

popt_10 = curve_fit(func_10, theta_10[:, 0], theta_10[:, 1])[0]
popt_30 = curve_fit(func_30, theta_30[:, 0], theta_30[:, 1])[0]

plt.plot(theta_10[:, 0], theta_10[:, 1])
plt.plot(theta_30[:, 0], theta_30[:, 1])

plt.plot(theta_10[:, 0], func_10(theta_10[:, 0], *popt_10))
plt.plot(theta_30[:, 0], func_30(theta_30[:, 0], *popt_30))

#%%
# V_0 = 3.06 eV
