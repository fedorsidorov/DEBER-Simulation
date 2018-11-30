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

#%% e_data
n_el = '3000'

e_data_source_dir = 'Wall/CUBES_C_exc_' + n_el + '/'

e_matrix = np.zeros((mv.n_Z_new, mv.n_XY, mv.n_x, mv.n_y, mv.n_z))

for k in range(mv.n_Z_new):
    
    print(k)
    
    for ij in range(mv.n_XY):
        
        e_data = np.load(e_data_source_dir + 'cube_' + str(k) + '_' + str(ij) + '.npy')
        
        for line in e_data:
            x, y, z = list(map(int, line))
            e_matrix[k, ij, x, y, z] += 1

#%
np.save('Wall/MATRIX_C_exc_' + n_el + '.npy', e_matrix)
#np.save('H_data/e_data_total_arr_ion.npy', e_data_total_arr)
