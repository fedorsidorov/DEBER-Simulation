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
## cube parameters
cube_size = 10.
cell_size = 2.

eps = 1e-3

## cubes number
n_total = mv.n_X * mv.n_Y * mv.n_Z
size = 1000

## arrays for future data
total_arr = np.zeros((n_total, size, 3)) * np.nan
total_list = [0] * n_total

## directories
source_dir = 'Sharma/DATA_Sharma/'
dest_dir = 'H_amount/CUBES_H_exc/'

fname_list = os.listdir(source_dir)

#%%
n_files = 0

for fname in fname_list:
    
    mf.upd_progress_bar(n_files, len(fname_list))
    
    if fname == '.DS_Store':
        continue
    
    DATA_Pn = np.load(source_dir + fname)
    
    ## apply statements to data
    statements = (np.all(DATA_Pn[:, (mv.e_x, mv.e_y, mv.e_z)] >= mv.xyz_min, axis=1),
                  np.all(DATA_Pn[:, (mv.e_x, mv.e_y, mv.e_z)] <= mv.xyz_max - eps, axis=1),
                  DATA_Pn[:, mv.coll_id] == mv.exc)
    
    inds = np.where(np.logical_and.reduce(statements))
    
    DATA_Pn_cut = DATA_Pn[inds]

    ## write DATA_Pn_cut to total_arr
    for line in DATA_Pn_cut:
        
        x, y, z = line[[mv.e_x, mv.e_y, mv.e_z]]
        xyz = line[[mv.e_x, mv.e_y, mv.e_z]]
        aid_Pn, dE_Pn = line[mv.atom_id], line[mv.e_dE]
        
        ## cut off extra dE
        if dE_Pn > 2e+4 or dE_Pn == 0:
            continue
        
        ## only H atoms are of interest
        if not aid_Pn == mv.H:
            continue
            
        ## x, y and z shifts in nm
        xyz_coords = xyz - mv.xyz_min
        
        ## get cube coordinates
        xyz_cube = np.floor(xyz_coords / cube_size).astype(int)
        
        ## get cube position in all_cubes_arr
        list_pos = np.dot(xyz_cube, (1, mv.n_X, mv.n_X * mv.n_Y)).astype(int)
        
        ## get cell coordinates in cube
        xyz_cell = np.floor((xyz_coords - xyz_cube * cube_size) / cell_size).astype(int)
        
        ## total
        total_arr[list_pos, total_list[list_pos], :] = xyz_cell
        total_list[list_pos] += 1
    
    n_files += 1

#%% write to files
for k in range(mv.n_Z):
    
    mf.upd_progress_bar(k, mv.n_Z)
    
    for ij in range(mv.n_XY):
    
        pos = k * (mv.n_X * mv.n_Y) + ij
        
        total_cube = mf.delete_nan_rows(total_arr[pos, :, :])
        
        np.save(dest_dir + 'cube_' + str(k) + '_' + str(ij) + '_total.npy', total_cube)
