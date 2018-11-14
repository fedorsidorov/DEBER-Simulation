#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
             
os.chdir(sim_path + 'mapping')

import sys
sys.path.append(sim_path + 'MODULES')

import importlib

import my_functions as mf
mf = importlib.reload(mf)

#%%
## cube parameters
cube_size = 10.
#cell_size = 0.5
#cell_size = 1.
cell_size = 2.

eps = 1e-3

## min and max coordinates
x_min, x_max = 0., 100.
y_min, y_max = 0., 100.
z_min, z_max = 0., 160.

## cubes numbers
n_x = int(np.round((x_max - x_min) / cube_size))
n_y = int(np.round((y_max - y_min) / cube_size))
n_z = int(np.round((z_max - z_min) / cube_size))

## cubes number
n_xy = n_x * n_y
n_total = n_x * n_y * n_z

## array for future data
total_arr = np.zeros((n_total, 50000, 3)) * np.nan
cubes_pos_hist = [0] * n_total

## directories
source_dir = 'Sharma/DATA_Sharma/'
dest_dir = 'Sharma/CUBES_C_exc/'

#%%
fname_list = os.listdir(source_dir)

n_files = 0

for fname in fname_list:
    
    mf.upd_progress_bar(n_files, len(fname_list))
    
    if fname == '.DS_Store':
        continue
    
    DATA_Pn = np.load(source_dir + fname)
    
    ## apple statements to data
    statements = (  DATA_Pn[:, 5] >= x_min, DATA_Pn[:, 5] <= x_max - eps \
                  , DATA_Pn[:, 6] >= y_min, DATA_Pn[:, 6] <= y_max - eps \
                  , DATA_Pn[:, 7] >= z_min, DATA_Pn[:, 7] <= z_max - eps \
                  ## the process is excitation
                  , DATA_Pn[:, 3] == 1
                  )
    inds = np.where(np.logical_and.reduce(statements))
    DATA_Pn_cut = DATA_Pn[inds]

    ## write DATA_Pn_cut to total_arr
    for line in DATA_Pn_cut:
        
        x, y, z = line[5:8]            
        aid_Pn, dE_Pn = line[2], line[8]
        
        if dE_Pn > 2e+4:
            continue
        
        ## 0 - H, 1 - C, 2 - 0, H8 C5 O2
        
        ## ... of 8 H atoms or ... of 5 C atoms are of interest
#        if (aid_Pn == 0 and mf.random() <= 0.125) or (aid_Pn == 1 and mf.random() <= 0.4):
            
        ## ... of 5 C atoms are of interest
        if aid_Pn == 1 and mf.random() <= 0.6:
            
            ## x, y and z shifts in nm
            x_coord = x - x_min
            y_coord = y - y_min
            z_coord = z - z_min
            
            ## get cube coordinates
            cube_x = np.floor(x_coord / cube_size)
            cube_y = np.floor(y_coord / cube_size)
            cube_z = np.floor(z_coord / cube_size)
            
            ## get cube position in all_cubes_arr
            cube_pos_1d = int(cube_x + cube_y * n_x + cube_z * n_x * n_y)
            
            ## get cell coordinates in cube
            cell_x = int(np.floor((x_coord - cube_x * cube_size) / cell_size))
            cell_y = int(np.floor((y_coord - cube_y * cube_size) / cell_size))
            cell_z = int(np.floor((z_coord - cube_z * cube_size) / cell_size))
            
            if cell_x < 0 or cell_y < 0 or cell_z < 0:
                print('!!!!')
            
            ## write cell coordinates into total_arr
            total_arr[cube_pos_1d, cubes_pos_hist[cube_pos_1d], :] =\
                cell_x, cell_y, cell_z
            cubes_pos_hist[cube_pos_1d] += 1
            
    n_files += 1

#%%
## write to files
for k in range(n_z):
    
    mf.upd_progress_bar(k, n_z)
    
    for ij in range(n_xy):
    
        pos = k * (n_x * n_y) + ij
        cube_arr = mf.delete_nan_rows(total_arr[pos, :, :])
        np.save(dest_dir + 'cube_' + str(k) + '_' + str(ij), cube_arr)
