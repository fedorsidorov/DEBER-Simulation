#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import copy
from itertools import product

import my_functions as mf
import my_variables as mv
import my_indexes as mi
import my_constants as mc
import my_mapping as mm

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mi = importlib.reload(mi)
mc = importlib.reload(mc)
mm = importlib.reload(mm)

os.chdir(mv.sim_path_MAC + 'maping')

#%%
e_matrix = np.load(mv.sim_path_MAC + 'MATRIXES/MATRIX_e_500_pC_cm_C.npy')
resist_matrix = np.load(mv.sim_path_MAC + 'MATRIXES/MATRIX_resist.npy')
chain_table = np.load(mv.sim_path_MAC + 'MATRIXES/TABLE_chains.npy')

p_scission = 0.5

n_scissions = 0

#%%
for x_ind, y_ind, z_ind in product(range(mm.resist_shape[0]),\
        range(mm.resist_shape[1]), range(mm.resist_shape[2])):
    
    if y_ind == z_ind == 0:
        mf.upd_progress_bar(x_ind, s_0)
    
    n_events = mm.get_n_events(e_matrix, x_ind, y_ind, z_ind)
    
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
            elif next_mon_type == 1:
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
        
        elif mon_type == mm.free_rad_mon:
            continue
        
        else:
            print('WTF', mon_type)

#%%
chain_test_inds = np.random.choice(mm.N_chains_total, 100, replace=False)
chain_test_table = chain_table[chain_test_inds]

resist_matrix_test = resist_matrix[:, 25, 30]

#%%
L_final = []
radical_matrix = np.zeros((len(chain_table), len(chain_table[0])))

for i, now_chain in enumerate(chain_table):
    
    mf.upd_progress_bar(i, N_chains_total)
    cnt = 0
    
    radical_matrix[i] = now_chain[:, mi.mon_type]
    
    for line in now_chain:
        
        if np.all(line == mc.uint16_max):
            break
        
        mon_type = line[-1]
                
        if mon_type in [mc.uint16_max, 9]:
            cnt == 1
        
        elif mon_type in [0, 10]:
            cnt += 1
        
        elif mon_type in [1, 11]:
            cnt += 1
            L_final.append(cnt)            
            cnt = 0

#%%
np.save('MATRIX_radicals.npy', radical_matrix)
np.save('final_L_2.5C_exc.npy', np.array(L_final))

#%%
radical_matrix_uint16 = np.array(radical_matrix, dtype=np.uint16)

#%%
L_final_arr = np.array(L_final)

#%%
plt.hist(np.log(L_final_arr * 100))

#%% get monomers
monomer_matrix = np.zeros((s_0, s_1, s_2))

for x_ind, y_ind, z_ind in product(range(s_0), range(s_1), range(s_2)):
    
    if y_ind == z_ind == 0:
        mf.upd_progress_bar(x_ind, s_0)
    
    now_cube = resist_matrix[x_ind, y_ind, z_ind]
    
    inds_2 = np.where(now_cube[:, mi.mon_type] == 2)[0]
    inds_12 = np.where(now_cube[:, mi.mon_type] == 12)[0]
    
    monomer_matrix[x_ind, y_ind, z_ind] += len(inds_2) + len(inds_12)

#%% drawing
plt.figure()
plt.semilogy(x_grid_2nm, np.sum(monomer_matrix[:, 25, :], axis = 1), label='500 pC/cm')
plt.xlabel('x, nm')
plt.ylabel('N monomers')
plt.title('Monomer coordinate distribution, 2 nm')
plt.legend()
plt.grid()
plt.show()
plt.savefig('LOG monomers 2nm.png', dpi=300)

#%%
monomer_matrix_xz = np.sum(monomer_matrix, axis=1)

plt.figure()
plt.semilogy(x_grid_2nm, np.sum(monomer_matrix_xz, axis = 1), label='500 pC/cm')
plt.xlabel('x, nm')
plt.ylabel('N monomers')
plt.title('Monomer coordinate distribution, 100 nm')
plt.legend()
plt.grid()
plt.show()
plt.savefig('LOG monomers 100nm.png', dpi=300)

#%% Sharma G-value
#N_el_dep = 6e-5 / 1.6e-19 * 1e-10
#E_dep =  N_el_dep * 25e+3
#G_value = n_scission / (E_dep / 100)
#print(G_value * 100)

#%%


#%% destroy some ester groups
#ester_part = 0.5
#
#add_CO = n_ester * ester_part * d_CO
#n_CO_final = n_CO + add_CO
#
#add_CO2 = n_ester * ester_part * d_CO2
#n_CO2_final = n_CO2 + add_CO2
#
#n_CH4_final = n_CH4 + (add_CO + add_CO2)
#n_ester_final = n_ester - (add_CO + add_CO2)
#
##%%
#print('n_CO', n_CO)
#print('n_CO2', n_CO2)
#print('n_CH4', n_CH4)
#print('n_ester', n_ester)
