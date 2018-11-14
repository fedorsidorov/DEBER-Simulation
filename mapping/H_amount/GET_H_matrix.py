#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import copy
import sys
import importlib

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'

sys.path.append(sim_path + 'MODULES')
import my_functions as mf
import my_variables as mv
mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(sim_path + 'mapping')

#%% cube parameters
cube_size = 10.
cell_size = 2.

## min and max coordinates
x_min, x_max = 0., 100.
xyz_min = np.array((0., 0., 0.))
xyz_max = np.array((100., 100., 160.))

## cubes parameters
n_cubes_x, n_cubes_y, n_cubes_z = (np.round((xyz_max - xyz_min) / cube_size)).astype(int)
n_cubes_xy = n_cubes_x * n_cubes_y
n_cell_x = n_cell_y = n_cell_z = int(cube_size / cell_size)

#%% e_data
e_data_source_dir = 'H_amount/CUBES_H_exc/'

H_matrix = np.zeros((n_cubes_z, n_cubes_xy, n_cell_x, n_cell_y, n_cell_z))

for k in range(n_cubes_z):
    
    print(k)
    
    for ij in range(n_cubes_xy):
        
        e_data = np.load(e_data_source_dir + 'cube_' + str(k) + '_' + str(ij) + '_total.npy')
        
        for line in e_data:
            x, y, z = list(map(int, line))
            H_matrix[k, ij, x, y, z] += 1

#%%
n_H = np.sum(H_matrix)


