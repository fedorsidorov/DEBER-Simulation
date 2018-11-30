#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv
import copy
from itertools import product
from scipy.integrate import odeint

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'mapping/radicals')

#%%
#n_rad_0 = 3740
#n_rad_0 = 37
n_rad_0 = 2721

N_a = 6.02e+23 # 1 / mole
V_cm3 = 4.8e-15 # cm^3
V_l = V_cm3 / 1000

R_0 = n_rad_0 / N_a / V_l

#%%
def model(y, t):
    
    k_t = 2850
#    k_t = 720 # litre mole-1 s-1
#    k_d = 68 # s-1
    
    dydt = - k_t * y**2 #- k_d * y
    
    return dydt


# time points
t = np.linspace(0, 3*3600, num=3600*3)

# solve ODE
y = odeint(model, R_0, t)

y_rad = y * N_a * V_l

# plot results
plt.plot(t, y_rad)
plt.xlabel('time, s')
plt.ylabel('[R]')
plt.grid()
#plt.plot(t, np.ones(len(t))*n_rad_0*0.6)
#plt.loglog(t, np.ones(len(t))*n_rad_0*0.01)
#plt.loglog(t, np.ones(len(t))*n_rad_0*0.01*0.6)
plt.show()

#%%
k_d = 68 # s-1
ind_1 = 1008
#ind_1 = 982

#k_d = 10
#ind_1 = 81

mon_conc = 0 # mole l-1

while ind_1 > 0:
    
    mon_conc += k_d * y[ind_1, 0]
    ind_1 -= 1

mon_n = mon_conc * N_a * V_l
mon_mol_cm3 = mon_conc / 1000

mon_mol_g = mon_mol_cm3 / 1.19


