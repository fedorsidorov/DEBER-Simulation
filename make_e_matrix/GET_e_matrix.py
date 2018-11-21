#%% Import
import numpy as np
import os
import importlib
import matplotlib.pyplot as plt

import my_functions as mf
import my_variables as mv
import e_matrix_functions as emf

mf = importlib.reload(mf)
mv = importlib.reload(mv)
emf = importlib.reload(emf)
os.chdir(mv.sim_path_MAC + 'make_e_matrix')

#%% calculate required files amount
n_files = emf.get_n_files_with_50nm_borders(500e-12, 100)

l_xyz = np.array((600, 100, 122))

space = 50

x_beg, y_beg, z_beg = (-l_xyz[0]/2, 0, 0)
xyz_beg = np.array((x_beg, y_beg, z_beg))
xyz_end = xyz_beg + l_xyz
x_end, y_end, z_end = xyz_end

step_2nm = 2

x_bins_2nm = np.arange(x_beg, x_end + 1, step_2nm)
y_bins_2nm = np.arange(y_beg, y_end + 1, step_2nm)
z_bins_2nm = np.arange(z_beg, z_end + 1, step_2nm)

bins_2nm = x_bins_2nm, y_bins_2nm, z_bins_2nm

#%%
n_mon_max = 400
e_matrix = np.zeros((len(x_bins_2nm)-1, len(y_bins_2nm)-1, len(z_bins_2nm)-1))

source_dir = mv.sim_path_MAC + 'e_DATA/DATA_Pn_20keV_122nm/'

n_files_used = 625

n_events = 0

for i in range(n_files_used):
    
    mf.upd_progress_bar(i, n_files_used)
    
    now_DATA_Pn = np.load(source_dir + 'DATA_Pn_' + str(i) + '.npy')
    
    n_rotations = 0
    
    for i in range(n_rotations + 1):
        
        if n_rotations > 0:
            emf.rotate_DATA(now_DATA_Pn)
        
        n_events += len(now_DATA_Pn)
        
        emf.shift_DATA(now_DATA_Pn, (x_beg, x_end), (-space + y_beg, y_end + space))
        
        e_matrix += np.histogramdd(now_DATA_Pn[:, 5:8], bins=bins_2nm)[0]
        
#%%
e_matrix_uint8 = np.array(e_matrix, dtype=np.uint8)
print('e_matrix size, Mb:', e_matrix_uint8.nbytes / 1024**2)
np.save('MATRIX_e_500_pC_cm.npy', e_matrix_uint8)  

