#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv
import copy
from itertools import product

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'mapping')

#%% total arrays
#e_matrix = np.load('Sharma/MATRIX_e_C_exc.npy')
e_matrix = np.load('Courtney/MATRIX_C_exc.npy')

part_matrix_0 = np.load('Sharma/MATRIX_particles.npy')
chain_inv_matrix_0 = np.load('Sharma/MATRIX_chains_inv.npy')

#%% Courtney
part_matrix = np.zeros((mv.n_Z_new, mv.n_XY, mv.n_x, mv.n_y, mv.n_z, mv.n_mon_max, 3))*np.nan
chain_inv_matrix = np.zeros((mv.N_chains*3, mv.chain_len_max, 7))*np.nan

for i in range(3):
    
    print(i)
    
    part_matrix[mv.n_Z * i : mv.n_Z * (i+1)] =\
        part_matrix_0[:, :, :, :, :, :mv.n_mon_max, :] 
    chain_inv_matrix[mv.N_chains * i : mv.N_chains * (i+1)] = chain_inv_matrix_0

#%%
CO_matrix = np.zeros(np.shape(e_matrix))
CO2_matrix = np.zeros(np.shape(e_matrix))
CH4_matrix = np.zeros(np.shape(e_matrix))

#%%
s0, s1, s2, s3, s4 = np.shape(e_matrix)

chain_lens = np.load('Sharma/initial_L.npy')

bond_scissions = [[] for i in range(mv.N_chains)]

## probabilities
p1 = 0.3
p2 = 0.5
p3 = 0.7

d1, d2, d3 = p1 - 0, p2 - p1, p3 - p2 

## scission ways
k2 = 4.6
k3 = 1.29

d_k2, d_k3 = (k2, k3) / np.sum((k2, k3))


for Z, XY, z, y, x in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
    
    if XY == z == y == x == 0:
        print(Z)
    
    n_events = e_matrix[Z, XY, x, y, z].astype(int)
    
    for i in range(n_events):
        
        ## only ... C atoms of 5 are of interest
        if mf.random() >= p3:
            continue
        
        ## current cell of working array 
        cell_part = part_matrix[Z, XY, z, y, x]
        
        ## indexes of all particles we have - monomers and ester groups
        part_inds = np.where(np.logical_not(np.isnan(cell_part[:, 0])))[0]
        
        if len(part_inds) == 0:
            continue
        
        now_part_ind = mf.choice(part_inds)
        now_part_line = cell_part[now_part_ind, :]
        
        
        
###############################################################################
# 0 - ester group #############################################################
###############################################################################
        if np.all(now_part_line == -1):
            
            ## remove ester group from the matrix
            cell_part[now_part_ind] = np.nan
            
            ## add one CO2_CH4 entry
            CH4_matrix[Z, XY, x, y, z] += 1
            
            easter_decay = mf.choice((mv.ester_CO, mv.ester_CO2), p=(d_k2, d_k3))
            
            if easter_decay == mv.ester_CO:
                CO_matrix[Z, XY, x, y, z] += 1
            else:
                CO2_matrix[Z, XY, x, y, z] += 1
            
            continue
    
    
    
        n_chain, n_mon, mon_type = list(map(int, now_part_line))
 


###############################################################################       
# 1, 2 - bonded monomer with or w/o ester group ###############################
############################################################################### 
        if mon_type == 0 or mon_type == 10:
            
            if mon_type == 0:
                probs = (d1, d2, d3)/np.sum((d1, d2, d3))
            else:
                probs = (0, 1, 0)
            
            process = mf.choice((mv.sci_ester, mv.sci_direct, mv.ester), p=probs)
            
            if process == mv.ester:
                
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -1] = 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
            else: ## deal with monomer types
                
                new_mon_kind = mf.choice((-1, 1))
                now_part_line[-1] += new_mon_kind
                chain_inv_matrix[n_chain, n_mon, -1] += new_mon_kind
                
                n_next_mon = n_mon + new_mon_kind
                next_mon_inv_line = chain_inv_matrix[n_chain, n_next_mon]
                
                Zn, XYn, xn, yn, zn = next_mon_inv_line[:5].astype(int)
                
                next_mon_pos = next_mon_inv_line[5]
                next_mon_type = next_mon_inv_line[6]
                
                ## if next monomer was at the end
                if next_mon_type == -1 or next_mon_type == -1 or\
                    next_mon_type == 9 or next_mon_type == 11:
                    
                    ## it becomes free
                    if next_mon_type == -1 or next_mon_type == 1: ## it has ester group
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 2
                    else:
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 12
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
                        
                        if next_mon_type == -1 or next_mon_type == 1:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 2
                        else:
                            chain_inv_matrix[n_chain, n_next_mon, -2] = 12
                        
                ## next monomer is not in the end
                else:
                    
                    next_mon_new_type = next_mon_type - new_mon_kind
                    
                    ## it becomes half-bonded
                    chain_inv_matrix[n_chain, n_next_mon, -1] = next_mon_new_type
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
        
                        part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] =\
                            next_mon_new_type
                        
            if process == mv.sci_ester: ## scission with ester group deatachment
                        
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -2] += 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
                
                
###############################################################################
# 3, 4 - half-bonded monomer with or w/o ester group ##########################
###############################################################################
        elif np.abs(mon_type) == 1 or np.abs(mon_type - 10) == 1:
            
            if np.abs(mon_type) == 1:
                probs = (d1, d2, d3)/np.sum((d1, d2, d3))
            else:
                probs = (0, 1, 0)
            
            process = mf.choice((mv.sci_ester, mv.sci_direct, mv.ester), p=probs)
            
            if process == mv.ester:
                
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -1] += 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
            else: ## deal with monomer types
                
                mon_kind = mon_type
                
                if np.abs(mon_type) == 1:
                    new_mon_type = 2
                else:
                    mon_kind -= 10
                    new_mon_type = 12
                
                now_part_line[-1] = new_mon_type
                chain_inv_matrix[n_chain, n_mon, -1] = new_mon_type
                
                n_next_mon = n_mon - mon_kind
                next_mon_inv_line = chain_inv_matrix[n_chain, n_next_mon]
                
                Zn, XYn, xn, yn, zn = next_mon_inv_line[:5].astype(int)
                
                next_mon_pos = next_mon_inv_line[5]
                next_mon_type = next_mon_inv_line[6]
                
                ## if next monomer was at the end
                if next_mon_type == -1 or next_mon_type == -1 or\
                    next_mon_type == 9 or next_mon_type == 11:
                    
                    ## it becomes free
                    if next_mon_type == -1 or next_mon_type == 1: ## it has ester group
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 2
                    else:
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 12
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
                        
                        if next_mon_type == -1 or next_mon_type == 1:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 2
                        else:
                            chain_inv_matrix[n_chain, n_next_mon, -2] = 12
                        
                ## next monomer is not in the end
                else:
                    
                    next_mon_new_type = next_mon_type - new_mon_type
                    
                    ## it becomes half-bonded
                    chain_inv_matrix[n_chain, n_next_mon, -1] = next_mon_new_type
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
        
                        part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] =\
                            next_mon_new_type
        
            if process == mv.sci_ester: ## scission with ester group deatachment
                        
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -2] += 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
        
        
        
###############################################################################
# 5 - bonded monomer with ester group #########################################
###############################################################################
        elif mon_type == 2:
            
            ## only ester group deatachment is possible
            
            now_part_line[-1] += 10
            chain_inv_matrix[n_chain, n_mon, -2] += 10
            
            ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
            cell_part[ester_ind, :] = -1
            
###############################################################################
# 5 - bonded monomer w/o ester group ##########################################
###############################################################################
        else:
            
            continue

#%% Analysis
n_ester_0 = 0

n_2 = 0
n_12 = 0

n_9 = 0
n_11 = 0

for Z, XY, z, y, x in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
    
    if XY == z == y == x == 0:
        print(Z)
    
    ## current cell of working array 
    cell_part = part_matrix[Z, XY, z, y, x]
    
    part_inds = np.where(np.logical_not(np.isnan(cell_part[:, 0])))
    
    for line in cell_part[part_inds]:
        
        if np.all(line == -1):
            
            n_ester_0 += 1
        
        if line[-1] == 2: n_2 += 1
        
        if line[-1] == 12: n_12 += 1
        
        if line[-1] == 9: n_9 += 1
        
        if line[-1] == 11: n_11 += 1

#%% Get L distribution
L_final = []

n = 0

for chain in chain_inv_matrix:
    
    mf.upd_progress_bar(n, 12709)
    n += 1
    cnt = 0
    
    for line in chain:
        
        if np.all(np.isnan(line)):
            break
        
        mon_type = line[-1]
                
        if mon_type == -1:
            cnt == 0
        
        elif mon_type == 0:
            cnt += 1
        
        elif mon_type == 1:
            L_final.append(cnt)            
            cnt = 0

#%%
np.save('Sharma/final_L_new.npy', np.array(L_final))

#%%
n_CO_0 = np.sum(CO_matrix)
n_CO2_0 = np.sum(CO2_matrix)
n_CH4_0 = np.sum(CH4_matrix)

#%% destroy some ester groups
ester_part = 0.5

add_CO = n_ester_0 * ester_part * d_k2
n_CO = n_CO_0 + add_CO

add_CO2 = n_ester_0 * ester_part * d_k3
n_CO2 = n_CO2_0 + add_CO2

n_CH4 = n_CH4_0 + (add_CO + add_CO2)
n_ester = n_ester_0 - (add_CO + add_CO2)

#%%
print('n_CO', n_CO)
print('n_CO2', n_CO2)
print('n_CH4', n_CH4)
print('n_ester', n_ester)
