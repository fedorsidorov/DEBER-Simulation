#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
#import copy
from itertools import product

import my_functions as mf
import my_variables as mv
import my_indexes as mi
import my_constants as mc

import scipy

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mi = importlib.reload(mi)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'mapping')

import my_mapping as mm
mm = importlib.reload(mm)

#%%
e_matrix = np.load('../MATRIXES/Harris/MATRIX_Harris_100uC_C_ion.npy')
e_matrix += np.load('../MATRIXES/Harris/MATRIX_Harris_100uC_C_exc.npy')
                
resist_matrix = np.load('../MATRIXES/Harris/MATRIX_resist_Harris.npy')
chain_table   = np.load('../MATRIXES/Harris/TABLE_chains_Harris.npy')
dE_matrix     = np.load('../MATRIXES/Harris/MATRIX_dE_Harris_100uC.npy')

N_chains_total  = len(chain_table)
N_mon_chain_max = len(chain_table[0])

resist_shape = np.shape(resist_matrix)[:3]

scission_matrix = np.zeros(np.shape(e_matrix))

sci_per_mol_matrix = np.zeros(N_chains_total)

#%%
#chain_sum_len_matrix_before, n_chains_matrix_before =\
#    mm.get_local_chain_len(resist_shape, N_mon_chain_max, chain_table, N_chains_total)
#
##%%
#np.save('Harris_chain_sum_len_matrix_before.npy', chain_sum_len_matrix_before)
#np.save('Harris_n_chains_matrix_before.npy', n_chains_matrix_before)

#%%
#p_scission = 0.5
p_scission = 1
n_scissions = 0

n_events_total = 0

#%%
for x_ind, y_ind, z_ind in product(range(resist_shape[0]),\
           range(resist_shape[1]), range(resist_shape[2])):
    
    if y_ind == z_ind == 0:
        mf.upd_progress_bar(x_ind, resist_shape[0])
    
    n_events = mm.get_n_events(e_matrix, x_ind, y_ind, z_ind)
    
    n_events_total += n_events
    
    for i in range(n_events):
        
        if mf.random() >= p_scission:
            continue
        
        resist_part_ind = mm.get_resist_part_ind(resist_matrix, x_ind, y_ind, z_ind)
        
        if resist_part_ind == -1:
            continue
        
        n_chain, n_mon, mon_type = mm.get_resist_part_line(resist_matrix,\
                                            x_ind, y_ind, z_ind, resist_part_ind)
        
        sci_per_mol_matrix[n_chain] += 1
        
############################################################################### 
        if mon_type == mm.mid_mon: ## bonded monomer ##########################
###############################################################################
            
            new_mon_type = mm.get_mon_type()
            
            mm.rewrite_mon_type(resist_matrix, chain_table,\
                                n_chain, n_mon, new_mon_type)
            
            n_next_mon = n_mon + new_mon_type - 1
            
            next_x_ind, next_y_ind, next_z_ind, _, next_mon_type =\
                mm.get_chain_table_line(chain_table, n_chain, n_next_mon)
            
            ## if next monomer was at the end
            if next_mon_type in [mm.beg_mon, mm.end_mon]:
                mm.rewrite_mon_type(resist_matrix, chain_table,\
                                 n_chain, n_next_mon, mm.free_mon)
            
            ## if next monomer is full bonded
            elif next_mon_type == mm.mid_mon:
                next_mon_new_type = next_mon_type - (new_mon_type - 1)
                mm.rewrite_mon_type(resist_matrix, chain_table,\
                                 n_chain, n_next_mon, next_mon_new_type)
            
            else:
                print('error 1, next_mon_type =', next_mon_type)
                print(x_ind, y_ind, z_ind)
            
            n_scissions += 1
            scission_matrix[x_ind, y_ind, z_ind] += 1

###############################################################################
        elif mon_type in [mm.beg_mon, mm.end_mon]: ## half-bonded monomer #####
###############################################################################
            
            new_mon_type = mm.free_mon
            
            mm.rewrite_mon_type(resist_matrix, chain_table,\
                                n_chain, n_mon, new_mon_type)
            
            n_next_mon = n_mon - (mon_type - 1) ## minus, Karl!
            
            next_x_ind, next_y_ind, next_z_ind, _, next_mon_type =\
                mm.get_chain_table_line(chain_table, n_chain, n_next_mon)
            
            ## if next monomer was at the end
            if next_mon_type in [mm.beg_mon, mm.end_mon]:
                mm.rewrite_mon_type(resist_matrix, chain_table,\
                                 n_chain, n_next_mon, mm.free_mon)
            
            ## if next monomer is full bonded
            elif next_mon_type == mm.mid_mon:
                next_mon_new_type = next_mon_type + (mon_type - 1)
                mm.rewrite_mon_type(resist_matrix, chain_table,\
                                 n_chain, n_next_mon, next_mon_new_type)
                
            else:
                print('error 2', next_mon_type)
            
            n_scissions += 1
            scission_matrix[x_ind, y_ind, z_ind] += 1
            
###############################################################################
        elif mon_type == mm.free_mon: ## free monomer with ester group ########
###############################################################################
            
            ## only ester group deatachment is possible
            mm.rewrite_mon_type(resist_matrix, chain_table,\
                             n_chain, n_mon, mm.free_rad_mon)
            
            n_scissions += 1
            scission_matrix[x_ind, y_ind, z_ind] += 1
        
        elif mon_type == mm.free_rad_mon:
            continue
        
        else:
            print('WTF', mon_type)

#%% G-value
E_dep = np.sum(dE_matrix)

G_value = n_scissions / (E_dep / 100)
#G_value = n_events_total / (E_dep / 100)

print(G_value)

#%%
L_final = []

for i, now_chain in enumerate(chain_table):
    
    mf.upd_progress_bar(i, N_chains_total)
    cnt = 0
    
    for line in now_chain:
        
        if np.all(line == mm.uint16_max):
            break
        
        mon_type = line[mm.mon_type_ind]
                
        if mon_type == 0:
            cnt == 1
        
        elif mon_type == 1:
            cnt += 1
        
        elif mon_type == 2:
            cnt += 1
            L_final.append(cnt)            
            cnt = 0

L_final_arr = np.array(L_final)

#%%
L_final_arr = np.load('Harris/L_final_2C_all_Harris.npy')

#%%
log_mw = np.log10(L_final_arr * 100)
plt.hist(log_mw, bins=20, cumulative=False, label='sample', rwidth=0.8,\
         normed=True, alpha=0.6)

data_B = np.loadtxt(mv.sim_path_MAC + 'make_chains/harris1973_B.dat')

x_B = data_B[:, 0]
y_B = data_B[:, 1]

x_B_log = np.log10(x_B)
X = np.linspace(x_B_log[0], x_B_log[-1], 200)
Y = mf.log_interp1d(x_B_log, y_B)(X)

#Y = scipy.signal.medfilt(Y, 5)

X_diff = X[:-1]
Y_diff = np.diff(Y)

#plt.plot(np.log10(x_B), y_B, label='model')
plt.plot(X_diff, Y_diff/np.max(Y_diff), label='model')

plt.title('Harris chain mass distribution after exposure, 2C ion')
plt.xlabel('log(m$_w$)')
plt.ylabel('probability')
plt.xlim((1.5, 5.5))
plt.ylim((0, 1))
plt.legend()
plt.grid()
plt.show()

#%%
chain_sum_len_matrix_C_1, n_chains_matrix_C_1 =\
    mm.get_local_chain_len(resist_shape, N_mon_chain_max, chain_table, N_chains_total)

#%%
np.save('chain_sum_len_matrix_C_1.npy', chain_sum_len_matrix_C_1)
np.save('n_chains_matrix_C_1.npy', n_chains_matrix_C_1)

#%%
sci_avg = np.average(sci_per_mol_matrix)

g = sci_avg * 1.19e-21 * 6.02e+23 / (np.sum(dE_matrix) / (100*100*500) * 9.5e+5)

print(g)
