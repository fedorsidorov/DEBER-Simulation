#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import copy
import time
import importlib

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'             

import sys
sys.path.append(sim_path + 'MODULES')
import my_functions as mf
mf = importlib.reload(mf)
import my_variables as mv
mv = importlib.reload(mv)

os.chdir(sim_path + 'mapping')

#%%
def get_x0y0(lx, ly, space):
    return mf.uniform(-space, lx + space), mf.uniform(-space, ly + space)

#%% cube parameters
cube_size = 10.
cell_size = 2.

eps = 1e-3

## cubes number
n_total = mv.n_X * mv.n_Y * mv.n_Z_new
size = 30000


lx = 100
ly = 100
space = 2


## arrays for future data
total_arr = np.zeros((n_total, size, 3)) * np.nan
total_list = [0] * n_total

## directories
#source_dir = 'Sharma/DATA_Sharma/'
source_dir = '../make_e_data/DATA_5keV_600nm/'
dest_dir = 'Courtney/CUBES_C_exc/'

fname_list = os.listdir(source_dir)

#%%
n_files = 0

n_H = 0
n_C = 0

for fname in fname_list:
    
    mf.upd_progress_bar(n_files, len(fname_list))
    
    if fname == '.DS_Store':
        continue
    
    DATA_Pn_raw = np.load(source_dir + fname)
    
    ## apply statements to data: C or H and excitation
    statements = (
            DATA_Pn_raw[:, mv.atom_id] != mv.O,
            DATA_Pn_raw[:, mv.coll_id] == mv.exc
            )
    
    inds = np.where(np.logical_and.reduce(statements))
    
    DATA_Pn = DATA_Pn_raw[inds]
    
    ## get more tracks!
    for i in range(5):
        
        ## copy DATA
        now_DATA_Pn = DATA_Pn.copy()
        
        ## make rotation
        now_DATA_Pn[:, 5:7] =\
            mf.add_xy_rotation(now_DATA_Pn[:, 5:7], 2*np.pi*mf.random())
        
        ## add shift
        ## primary tracks number
        n_tr_prim = int(now_DATA_Pn[np.where(np.isnan(now_DATA_Pn[:, 1]))][-1, 0] + 1)
        
        for track_num in range(n_tr_prim):
        
            ## in case of only elastic events in PMMA
            if len(np.where(now_DATA_Pn[:, 0] == track_num)[0]) == 0:
                continue
            
            ## in normal case
            else:
                x0, y0 = get_x0y0(lx, ly, space)            
                now_DATA_Pn = mf.add_xy_shift(now_DATA_Pn, track_num, x0, y0)
        
        ## apply statements to data: x, y, z limits
        statements = (np.all(now_DATA_Pn[:, (mv.e_x, mv.e_y, mv.e_z)] >= mv.xyz_min, axis=1),
                      np.all(now_DATA_Pn[:, (mv.e_x, mv.e_y, mv.e_z)] <= mv.xyz_max_new - eps,\
                                         axis=1))
        
        inds = np.where(np.logical_and.reduce(statements))
        
        inds_H = np.where(now_DATA_Pn[:, mv.atom_id] == mv.H)[0]
        inds_C = np.where(now_DATA_Pn[:, mv.atom_id] == mv.C)[0]
        
        n_H += len(inds_H)
        n_C += len(inds_C)
        
        continue
        
        now_DATA_Pn_cut = now_DATA_Pn[inds]
        
        ## write DATA_Pn_cut to total_arr
        for line in now_DATA_Pn_cut:
            
            xyz = line[[mv.e_x, mv.e_y, mv.e_z]]
            
            aid_Pn, dE_Pn = line[mv.atom_id], line[mv.e_dE]
            
            ## cut off extra dE
            if dE_Pn > 2e+4:
                continue
                    
            ## x, y and z shifts in nm
            xyz_coords = xyz - mv.xyz_min
            
            ## get cube coordinates
            xyz_cube = np.floor(xyz_coords / cube_size).astype(int)
            
            ## get cube position in all_cubes_arr
            list_pos = np.dot(xyz_cube, (1, mv.n_X, mv.n_X * mv.n_Y)).astype(int)
            
            ## get cell coordinates in cube
            xyz_cell = np.floor((xyz_coords - xyz_cube * cube_size) / cell_size).astype(int)
            
            ## write to total array
            total_arr[list_pos, total_list[list_pos], :] = xyz_cell
            total_list[list_pos] += 1
            
    n_files += 1

#%% write to files
for k in range(mv.n_Z_new):
    
    mf.upd_progress_bar(k, mv.n_Z_new)
    
    for ij in range(mv.n_XY):
        
        pos = k * (mv.n_X * mv.n_Y) + ij
        
        total_cube = mf.delete_nan_rows(total_arr[pos, :, :])
        
        np.save(dest_dir + 'cube_' + str(k) + '_' + str(ij) + '.npy', total_cube)
