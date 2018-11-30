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

#%% total arrays - 15 s
#e_matrix = np.load('Sharma/MATRIX_e_C_exc.npy')
e_matrix = np.load('Courtney/MATRIX_C_exc.npy')

part_matrix_0 = np.load('Sharma/MATRIX_particles.npy')
chain_inv_matrix_0 = np.load('Sharma/MATRIX_chains_inv.npy')

#% load part_matrix - 30 s
part_matrix = np.zeros((mv.n_Z_new, mv.n_XY, mv.n_x, mv.n_y, mv.n_z, 500, 3))*np.nan

for i in range(3):
    
    print(i)
    part_matrix[mv.n_Z * i : mv.n_Z * (i+1)] =\
        part_matrix_0

#% make correction to part_matrix - 15 s
for i in range(1, 3):
    
    print(i)
    
    for i0 in range(mv.n_Z * i, mv.n_Z*(i + 1)):
        
        for i1, i2, i3, i4 in product(range(mv.n_XY), range(mv.n_x),\
                                      range(mv.n_y), range(mv.n_z)):
    
            part_matrix[i0, i1, i2, i3, i4] += mv.n_chains * i, 0, 0

#% load chain_inv_matrix - 90 s
chain_inv_matrix = np.zeros((mv.n_chains*3, mv.chain_len_max, 7))*np.nan

for i in range(3):
    
    print(i)
    chain_inv_matrix[mv.n_chains * i : mv.n_chains * (i+1)] = chain_inv_matrix_0

#% make correction to chain_inv_matrix - 45 s
for i in range(1, 3):
    
    print(i)
    
    for i0 in range(mv.n_chains * i, mv.n_chains*(i + 1)):
        chain_inv_matrix[i0] += mv.n_Z * i, 0, 0, 0, 0, 0, 0

#%
#CO_matrix = np.zeros(np.shape(e_matrix))
#CO2_matrix = np.zeros(np.shape(e_matrix))
#CH4_matrix = np.zeros(np.shape(e_matrix))

n_CO_0 = 0
n_CO2_0 = 0
n_CH4_0 = 0

#%% 16 min
s0, s1, s2, s3, s4 = np.shape(e_matrix)

chain_lens = np.load('Sharma/initial_L.npy')

## probabilities
p1 = 0.3 ## ester group detachment with scissions
p2 = 0.5 ## sure lol
p3 = 0.7 ## ester group detachment w/o scissions

d1, d2, d3 = p1 - 0, p2 - p1, p3 - p2 

## scission ways
k_CO = 25.3
k_CO2 = 13

d_CO, d_CO2 = (k_CO, k_CO2) / np.sum((k_CO, k_CO2))

n_ester_0 = 0

for Z, XY, z, y, x in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
#for Z, XY, z, y, x in product(range(0, 1), range(0, 1), range(0, 1), range(0, 1), range(0, 1)):
    
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
#        part_inds = np.where(np.logical_not(np.isnan(part_matrix[Z, XY, z, y, x, :, 0])))[0]
        
        if len(part_inds) == 0:
            continue
        
        now_part_ind = mf.choice(part_inds)
        now_part_line = cell_part[now_part_ind, :]
#        now_part_line = part_matrix[Z, XY, z, y, x, now_part_ind, :]
        
        
###############################################################################
# 0 - ester group #############################################################
###############################################################################
        if np.all(now_part_line == -1):
            
            ## remove ester group from the matrix
            cell_part[now_part_ind] = np.nan
#            part_matrix[Z, XY, z, y, x, now_part_ind, :] = np.nan
            
            ## add one CO2_CH4 entry
#            CH4_matrix[Z, XY, x, y, z] += 1
            n_CH4_0 += 1
            
            easter_decay = mf.choice((mv.ester_CO, mv.ester_CO2), p=(d_CO, d_CO2))
            
            if easter_decay == mv.ester_CO:
#                CO_matrix[Z, XY, x, y, z] += 1
                n_CO_0 += 1
            else:
#                CO2_matrix[Z, XY, x, y, z] += 1
                n_CO2_0 += 1
            
            n_ester_0 -= 1
            
            continue
    
    
    
        n_chain, n_mon, mon_type = list(map(int, now_part_line))
 


###############################################################################       
# 1, 2 - bonded monomer with or w/o ester group ###############################
############################################################################### 
        if mon_type in [0, 10]:
            
            if mon_type == 0:
                probs = (d1, d2, d3)/np.sum((d1, d2, d3))
            else:
                probs = (0, 1, 0)
            
            process = mf.choice((mv.sci_ester, mv.sci_direct, mv.ester), p=probs)
            
            if process == mv.ester:
                
                if now_part_line[-1] != 0:
                    print('error 0')
                    print(process)
                    print(now_part_line)
                
                now_part_line[-1] = 10
                
                if chain_inv_matrix[n_chain, n_mon, -1] != 0:
                    print('error 1')
                    print(process)
                    print(now_part_line)
                
                chain_inv_matrix[n_chain, n_mon, -1] = 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
#                ester_ind =\
#                    np.where(np.isnan(part_matrix[Z, XY, z, y, x, :, 0]))[0][0]
                cell_part[ester_ind, :] = -1
#                part_matrix[Z, XY, z, y, x, now_part_ind, :] = -1
                
                n_ester_0 += 1
                
            else: ## deal with monomer types
                
                ## new kind of monomer
                new_mon_kind = mf.choice([-1, 1])
                
                now_part_line[-1] += new_mon_kind
                chain_inv_matrix[n_chain, n_mon, -1] += new_mon_kind
                
                n_next_mon = n_mon + new_mon_kind
                next_mon_inv_line = chain_inv_matrix[n_chain, n_next_mon]
                
                Zn, XYn, xn, yn, zn = next_mon_inv_line[:5].astype(int)
                
                next_mon_pos = next_mon_inv_line[5]
                next_mon_type = next_mon_inv_line[6]
                
                if next_mon_type not in [-1, 0, 1, 9, 10, 11]:
                    print('error 2')
                    print('process:', process)
                    print('now part line:', now_part_line)
                    print('next_mon_type:', next_mon_type)
                    print('n_mon', n_mon)
                    print('n_next_mon', n_next_mon)
                
                ## if next monomer was at the end
                if next_mon_type in [-1, 1, 9, 11]:
                        
                    ## it becomes free
                    if next_mon_type in [-1, 1]: ## it has ester group
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 2
                    else:
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 12
                    
                    if not np.isnan(next_mon_pos): ## write to part matrix if needed
                        
                        if next_mon_type in [-1, 1]:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 2
                        else:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 12
                        
                ## next monomer is not in the end
                else:
                    
                    next_mon_new_type = next_mon_type - new_mon_kind
                    
                    ## it becomes half-bonded
                    chain_inv_matrix[n_chain, n_next_mon, -1] = next_mon_new_type
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
        
                        part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] =\
                            next_mon_new_type
                        
            if process == mv.sci_ester: ## scission with ester group deatachment
                
                now_part_line[-1] = 10
                chain_inv_matrix[n_chain, n_mon, -1] = 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
#                ester_ind =\
#                    np.where(np.isnan(part_matrix[Z, XY, z, y, x, :, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
                n_ester_0 += 1
                
#%%        
        continue
        
        print('here')
        
###############################################################################
# 3, 4 - half-bonded monomer with or w/o ester group ##########################
###############################################################################
        if mon_type in [-1, 1, 9, 11]:
            
            if mon_type in [-1, 1]:
                probs = (d1, d2, d3)/np.sum((d1, d2, d3))
            else:
                probs = (0, 1, 0)
            
            process = mf.choice((mv.sci_ester, mv.sci_direct, mv.ester), p=probs)
            
            if process == mv.ester:
                
                if now_part_line[-1] > 1:
                    print('error 3')
                    print(process)
                    print(now_part_line)
                
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -1] += 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
                n_ester_0 += 1
                
            else: ## deal with monomer types
                
                mon_kind = mon_type
                
                if mon_type in [-1, 1]:
                    new_mon_type = 2
                else:
                    mon_kind -= 10
                    new_mon_type = 12
                    
                if mon_kind not in [-1, 1]:
                    print('error 4')
                    print(process)
                    print(now_part_line)
                
                now_part_line[-1] = new_mon_type
                chain_inv_matrix[n_chain, n_mon, -1] = new_mon_type
                
                n_next_mon = n_mon - mon_kind
                next_mon_inv_line = chain_inv_matrix[n_chain, n_next_mon]
                
                Zn, XYn, xn, yn, zn = next_mon_inv_line[:5].astype(int)
                
                next_mon_pos = next_mon_inv_line[-2]
                next_mon_type = next_mon_inv_line[-1]
                
                if next_mon_type not in [-1, 0, 1, 9, 10, 11]:
                    print('error 4.5')
                    print(process)
                    print(now_part_line)
                
                ## next monomer is in the end
                if next_mon_type in [-1, 1, 9, 11]:
                        
                    ## it becomes free
                    if next_mon_type in [-1, 1]: ## it has ester group
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 2
                    else:
                        chain_inv_matrix[n_chain, n_next_mon, -1] = 12
                    
                    if not np.isnan(next_mon_pos): ## write to part matrix if needed
                        
                        if next_mon_type in [-1, 1]:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 2
                        else:
                            part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] = 12
                        
                ## next monomer is not in the end
                else:
                    
                    if next_mon_type not in [0, 10]:
                        print('error 4.5')
                        print(process)
                        print(now_part_line)
                    
                    next_mon_new_type = next_mon_type - mon_kind
                    
                    if next_mon_new_type not in [-1, 1, 9, 11]:
                        print('error 5')
                        print(process)
                        print(now_part_line)
                    
                    ## it becomes half-bonded
                    chain_inv_matrix[n_chain, n_next_mon, -1] = next_mon_new_type
                    
                    if not np.isnan(next_mon_pos): ## write  to part matrix if needed
        
                        part_matrix[Zn, XYn, xn, yn, zn, int(next_mon_pos), -1] =\
                            next_mon_new_type
        
            if process == mv.sci_ester: ## scission with ester group deatachment
                
                if now_part_line[-1] != 2:
                    print('error 6')
                    print(process)
                    print(now_part_line)
                
                now_part_line[-1] += 10
                chain_inv_matrix[n_chain, n_mon, -2] += 10
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
                cell_part[ester_ind, :] = -1
                
                n_ester_0 += 1
        
###############################################################################
# 5 - free monomer with ester group ###########################################
###############################################################################
        elif mon_type == 2:
            
            ## only ester group deatachment is possible
            
            now_part_line[-1] = 12
            chain_inv_matrix[n_chain, n_mon, -1] = 12
            
            ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0]
            cell_part[ester_ind, :] = -1
            
            n_ester_0 += 1
            
###############################################################################
# 5 - free monomer w/o ester group ############################################
###############################################################################
        else:
            
            if mon_type != 12:
                print('error 7')
                print(process)
                print(now_part_line)
            
            continue

#%% Analysis
n_ester_0 = 0
#
#n_2 = 0
#n_12 = 0
#
#n_9 = 0
#n_11 = 0
#
#for Z, XY, z, y, x in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
#    
#    if XY == z == y == x == 0:
#        print(Z)
#    
#    ## current cell of working array 
#    cell_part = part_matrix[Z, XY, z, y, x]
#    
#    part_inds = np.where(np.logical_not(np.isnan(cell_part[:, 0])))
#    
#    for line in cell_part[part_inds]:
#        
#        if np.all(line == -1):
#            
#            n_ester_0 += 1
#        
#        if line[-1] == 2: n_2 += 1
#        
#        if line[-1] == 12: n_12 += 1
#        
#        if line[-1] == 9: n_9 += 1
#        
#        if line[-1] == 11: n_11 += 1

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
#n_CO_0 = np.sum(CO_matrix)
#n_CO2_0 = np.sum(CO2_matrix)
#n_CH4_0 = np.sum(CH4_matrix)

#%% destroy some ester groups
ester_part = 0.5

add_CO = n_ester_0 * ester_part * d_CO
n_CO = n_CO_0 + add_CO

add_CO2 = n_ester_0 * ester_part * d_CO2
n_CO2 = n_CO2_0 + add_CO2

n_CH4 = n_CH4_0 + (add_CO + add_CO2)
n_ester = n_ester_0 - (add_CO + add_CO2)

#%%
print('n_CO', n_CO)
print('n_CO2', n_CO2)
print('n_CH4', n_CH4)
print('n_ester', n_ester)
