#%% Import
import numpy as np
import importlib
#import my_arrays as ma
import my_functions_no2 as mf
mf = importlib.reload(mf)
#import my_variables as mv
#import matplotlib.pyplot as plt

#%%
n_files = 500
n_tracks = 10
d_PMMA = 1e+5
D = 0
beam_pos = 0
num = 150

while num < n_files:
    
    DATA = mf.get_DATA(1e+6, beam_pos, D, d_PMMA, n_tracks)
    
    fname = 'DATA_PMMA_Si/DATA_PMMA_Si_' + str(num) + '.npy'
    np.save(fname, DATA)
    print('file ' + fname + ' is saved')

    num += 1

#%%
mf.plot_DATA(DATA, 1)
