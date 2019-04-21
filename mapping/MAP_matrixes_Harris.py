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

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mi = importlib.reload(mi)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'mapping')

import my_mapping as mm
mm = importlib.reload(mm)

#%%
e_matrix =      np.load('../MATRIXES/Harris/MATRIX_Harris_100uC_C_exc.npy')
resist_matrix = np.load('../MATRIXES/Harris/MATRIX_resist_Harris.npy')
chain_table =   np.load('../MATRIXES/Harris/TABLE_chains_Harris.npy')

dE_matrix =     np.load('../MATRIXES/Harris/MATRIX_dE_Harris_100uC.npy')

N_chains_total = len(chain_table)
N_mon_chain_max = len(chain_table[0])

resist_shape = np.shape(resist_matrix)[:3]

#%%
#chain_sum_len_matrix_before, n_chains_matrix_before =\
#    mm.get_local_chain_len(resist_shape, N_mon_chain_max, chain_table, N_chains_total)
#
##%%
#np.save('Harris_chain_sum_len_matrix_before.npy', chain_sum_len_matrix_before)
#np.save('Harris_n_chains_matrix_before.npy', n_chains_matrix_before)

#%%
p_scission = 0.5
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
            
###############################################################################
        elif mon_type == mm.free_mon: ## free monomer with ester group ########
###############################################################################
            
            ## only ester group deatachment is possible
            mm.rewrite_mon_type(resist_matrix, chain_table,\
                             n_chain, n_mon, mm.free_rad_mon)
            
            n_scissions += 1
        
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


L_final = mm.get_L_final(chain_table)

#%%
chain_sum_len_matrix_C_1, n_chains_matrix_C_1 =\
    mm.get_local_chain_len(resist_shape, N_mon_chain_max, chain_table, N_chains_total)

#%%
np.save('chain_sum_len_matrix_C_1.npy', chain_sum_len_matrix_C_1)
np.save('n_chains_matrix_C_1.npy', n_chains_matrix_C_1)
