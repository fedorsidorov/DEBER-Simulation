#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'test_chain_size')

#%%
match = 0
mismatch = 0

R = 3

for i in range(1000):

    now_arr = np.load(mv.sim_path_MAC + 'test_chain_size/160_mon/chain_' + str(i) + '.npy')
    
    xc, yc, zc = np.average(now_arr, axis=0)
    
    for mon in now_arr:
        
        if (mon[0] - xc)**2 + (mon[1] - yc)**2 + (mon[2] - zc)**2 <= R**2:
            match += 1
        
        else:
            mismatch += 1

#%%
mon_mass = (12*5 + 16*2 + 1*8) / 6.02e+23

n_mon_4nm = 4/3 * np.pi * (2e-7)**3 * 1.19 / mon_mass