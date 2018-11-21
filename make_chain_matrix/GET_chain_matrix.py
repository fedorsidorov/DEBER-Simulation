#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'make_chain_matrix')

#%%
source_dir = mv.sim_path_MAC + 'CHAINS/950K_122nm/comb_600x100x122_center/'

N_0 = 63306
max_len = 9780

l_xyz = np.array((600, 100, 122))

x_beg, y_beg, z_beg = (-l_xyz[0]/2, 0, 0)
xyz_beg = np.array((x_beg, y_beg, z_beg))
xyz_end = xyz_beg + l_xyz
x_end, y_end, z_end = xyz_end

step_2nm = 2

x_bins_2nm = np.arange(x_beg, x_end + 1, step_2nm)
y_bins_2nm = np.arange(y_beg, y_end + 1, step_2nm)
z_bins_2nm = np.arange(z_beg, z_end + 1, step_2nm)

n_mon_max = 400

pos_matrix = np.zeros(l_xyz, dtype=np.uint16)
chain_matrix = - np.ones((*l_xyz, n_mon_max, 3), dtype=np.uint16)

#%%
for chain_num in range(N_0):
    
    mf.upd_progress_bar(chain_num, N_0)
    
    now_chain = np.load(source_dir + 'chain_shift_' + str(chain_num) + '.npy')
    
    for mon_pos, mon_line in enumerate(now_chain):
        
        if not (np.all(mon_line >= xyz_beg) and np.all(mon_line <= xyz_end)):
            continue
        
        now_x, now_y, now_z = mon_line
        
        x_ind = mf.get_closest_el_ind(x_bins_2nm, now_x)
        y_ind = mf.get_closest_el_ind(y_bins_2nm, now_y)
        z_ind = mf.get_closest_el_ind(z_bins_2nm, now_z)
        
        if mon_pos == 0:
            mon_type = -1
        elif mon_pos == len(now_chain) - 1:
            mon_type = 1
        else:
            mon_type = 0
        
        chain_matrix[x_ind, y_ind, z_ind, pos_matrix[x_ind, y_ind, z_ind]] =\
            chain_num, mon_pos, mon_type
        
        pos_matrix[x_ind, y_ind, z_ind] += 1

#%%
print('chain_matrix size, Gb:', chain_matrix.nbytes / 1024**3)
np.save('MATRIX_chain.npy', chain_matrix)        

#%%
cim = np.load('MATRIX_chain_inv.npy')

#%%
pm = pos_matrix[0]
cm = chain_matrix[0, 0, 3]




