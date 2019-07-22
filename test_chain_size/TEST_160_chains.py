#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
#import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'test_chain_size')

#%%
match = 0
mismatch = 0

R = 2

for i in range(10000):

    chain_arr = np.load(mv.sim_path_MAC + 'CHAINS/160_mon/chain_'+ str(i) + '.npy')[:160]
    
    xyz_c = np.average(chain_arr, axis=0)
    
    chain_arr_shift = chain_arr - xyz_c
    
    for mon in chain_arr_shift:
        
        if np.linalg.norm(mon) <= R:
            match += 1
        else:
            mismatch += 1
    
    
    
    
    
    
