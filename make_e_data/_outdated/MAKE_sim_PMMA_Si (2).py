#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'make_e_data')

#%%
## Usual
n_files = 1000
n_tracks = 100

## Bruk
d_PMMA = 80
E0 = 20e+3

## Sharma
#d_PMMA = 160
#E0 = 25e+3

## Cortney
#d_PMMA = 600
#E0 = 5e+3

## Aktary, dose = 100 uC/cm^2, A = 20x20 nm^2
## n_e = 100e-6 * (20e-7)**2 / 1.6e-19 = 2500
#d_PMMA = 100
#E0 = 10e+3

## CASINO
#d_PMMA = 1000
#E0 = 20e+3

D = 0
num = 24

while num < n_files:
#while num < 1:
    
    DATA = mf.get_DATA(E0, D, d_PMMA, n_tracks)
    
    DATA_P = DATA[np.where(DATA[:, mv.atom_id] != mv.Si)[0], :]
    DATA_Pn = DATA_P[np.where(DATA_P[:, mv.coll_id] > mv.elastic)[0], :]
    
    fname = 'DATA_Pn_20keV_80nm/DATA_Pn_' + str(num) + '.npy'
    np.save(fname, DATA_Pn)
    
    print('file ' + fname + ' is saved')

    num += 1

#%%
DATA = np.load('DATA_20keV_80nm/DATA_0.npy')

#%%
mf.plot_DATA(DATA, d_PMMA)
