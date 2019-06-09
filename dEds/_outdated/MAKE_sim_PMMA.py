#%% Import
import numpy as np
#import importlib
#import my_simulation_3 as ms
#ms = importlib.reload(ms)
#import my_arrays as ma
import sys
sys.path.append('../MODULES')
import my_functions as mf
#import my_variables as mv
#import matplotlib.pyplot as plt

#%%
n_files = 100
n_tracks = 50
d_PMMA = 1e0
E0 = 25e+3
D = 0
num = 2

while num < n_files:
    
    DATA = mf.get_DATA(E0, D, d_PMMA, n_tracks)
    
#    DATA_P = DATA[np.where(DATA[:, 2] != 3)[0], :]
#    DATA_Pn = DATA_P[np.where(DATA_P[:, 3] > 0)[0], :]
    
    fname = 'DATA/DATA_PMMA_' + str(num) + '.npy'
    np.save(fname, DATA)
    print('file ' + fname + ' is saved')

    num += 1

#%%
DATA = np.load('DATA/DATA_PMMA_0.npy')

mf.plot_DATA(DATA)
