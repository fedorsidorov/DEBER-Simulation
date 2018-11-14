#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv
import copy
from itertools import product

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'mapping')

#%%
n_el_arr = 100, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000
n_radicals_arr = []

for ind, n_el in enumerate(n_el_arr[3:4]):
    
    mf.upd_progress_bar(ind, len(n_el_arr))
    
    radical_matrix = np.load('Wall/radical_matrix_' + str(n_el) + '.npy')
    
    n_radicals = 0
    
    for line in radical_matrix:

        n_radicals_now = len(np.where(line == 1)[0]) + len(np.where(line == 11)[0])
        n_radicals += n_radicals_now -1
    
    n_radicals_arr.append(n_radicals)

#%%
plt.plot(n_el_arr, n_radicals_arr, 'ro')
plt.xlabel('n_electrons')
plt.ylabel('n_radicals')

#%%
N_a = 6.02e+23
V = 4.8e-15

t_arr = np.array(n_el_arr) / 500
n_rad_moles_cm3 = np.array(n_radicals_arr) / N_a / V

plt.figure()

plt.plot(t_arr, n_rad_moles_cm3, 'o-', label='my radicals')

#%
paper_radicals = np.loadtxt('radicals/Onishi_radicals.csv')
xdata = paper_radicals[:, 0]
ydata = paper_radicals[:, 1]*1e-7

plt.plot(xdata, ydata, 'o-', label='paper')
plt.xlabel('time, hours')
plt.ylabel('n_radicals')
plt.legend()

ax = plt.gca()
#ax.ticklabel_format(style='sci')
ax.yaxis.get_major_formatter().set_powerlimits((0, 1))

plt.show()

