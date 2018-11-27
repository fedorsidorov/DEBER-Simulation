#%% Import
import numpy as np
import os
import importlib

import my_functions as mf
import my_variables as mv
import my_constants as mc
import my_indexes as mi

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mc = importlib.reload(mc)
mi = importlib.reload(mi)

os.chdir(mv.sim_path_MAC + 'make_chain_matrix')

#%%
source_dir = mv.sim_path_MAC + 'CHAINS/950K_122nm/comb_600x100x122_center/'

N_0 = 63306
max_len = 9780
#n_mon_max = 404
n_mon_max = 810

l_xyz = np.array((600, 100, 122))

x_beg, y_beg, z_beg = (-l_xyz[0]/2, 0, 0)
xyz_beg = np.array((x_beg, y_beg, z_beg))
xyz_end = xyz_beg + l_xyz
x_end, y_end, z_end = xyz_end

step_2nm = 2

x_bins_2nm = np.arange(x_beg, x_end + 1, step_2nm)
y_bins_2nm = np.arange(y_beg, y_end + 1, step_2nm)
z_bins_2nm = np.arange(z_beg, z_end + 1, step_2nm)

x_grid_2nm = (x_bins_2nm[:-1] + x_bins_2nm[1:]) / 2
y_grid_2nm = (y_bins_2nm[:-1] + y_bins_2nm[1:]) / 2
z_grid_2nm = (z_bins_2nm[:-1] + z_bins_2nm[1:]) / 2

#%%
shape_arr = len(x_grid_2nm), len(y_grid_2nm), len(z_grid_2nm)

pos_matrix = np.zeros(shape_arr, dtype=np.uint16)
resist_matrix = - np.ones((*shape_arr, n_mon_max, 3), dtype=np.uint16)

chain_table = - np.ones((N_0, max_len, 5), dtype=np.uint16)

#%%
for chain_num in range(N_0):
#for chain_num in range(10):
    
    mf.upd_progress_bar(chain_num, N_0)
    
    now_chain = np.load(source_dir + 'chain_shift_' + str(chain_num) + '.npy')
    
    for n_mon, mon_line in enumerate(now_chain):
        
        if n_mon == 0:
            mon_type = mc.uint16_max
        elif n_mon == len(now_chain) - 1:
            mon_type = 1
        else:
            mon_type = 0
        
        if not (np.all(mon_line >= xyz_beg) and np.all(mon_line <= xyz_end)):
            chain_table[chain_num, n_mon, mi.mon_type] = mon_type
            continue
        
        now_x, now_y, now_z = mon_line
        
        x_ind = mf.get_closest_el_ind(x_grid_2nm, now_x)
        y_ind = mf.get_closest_el_ind(y_grid_2nm, now_y)
        z_ind = mf.get_closest_el_ind(z_grid_2nm, now_z)
        
        mon_line_pos = pos_matrix[x_ind, y_ind, z_ind]
        
        resist_matrix[x_ind, y_ind, z_ind, mon_line_pos] = chain_num, n_mon, mon_type
        
        chain_table[chain_num, n_mon] = x_ind, y_ind, z_ind, mon_line_pos, mon_type
        
        pos_matrix[x_ind, y_ind, z_ind] += 1
        
#%%        
print('resist_matrix size, Gb:', resist_matrix.nbytes / 1024**3)
np.save('MATRIX_resist.npy', resist_matrix)

#%%
print('chain_table size, Gb:', chain_table.nbytes / 1024**3)
np.save('TABLE_chains.npy', chain_table)
