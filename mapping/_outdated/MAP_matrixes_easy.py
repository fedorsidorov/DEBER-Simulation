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
e_matrix = np.load('Sharma/MATRIX_e_C_exc.npy')
part_matrix = np.load('Sharma/MATRIX_easy_mon.npy')

ester_matrix = np.zeros(np.shape(e_matrix))

CO_matrix = np.zeros(np.shape(e_matrix))
CO2_matrix = np.zeros(np.shape(e_matrix))
CH4_matrix = np.zeros(np.shape(e_matrix))

#%%
s0, s1, s2, s3, s4 = np.shape(e_matrix)

## probabilities
p1 = 0.3
p2 = 0.5
p3 = 0.7

d1, d2, d3 = p1 - 0, p2 - p1, p3 - p2

k2 = 4.6
k3 = 1.29

d_k2, d_k3 = (k2, k3) / np.sum((k2, k3))

for Z, XY, x, y, z in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
    
    if XY == z == y == x == 0:
        print(Z)
    
    n_events = e_matrix[Z, XY, x, y, z].astype(int)
    
    for i in range(n_events):
        
        ## only ... C atoms of 5 are of interest
        if mf.random() >= p3:
            continue
        
        ## current cell of working array 
        cell_part = part_matrix[Z, XY, x, y, z]
        
        ## indexes of all particles we have - monomers and ester groups
        part_inds = np.where(np.logical_not(np.isnan(cell_part[:, 0])))[0].astype(int)
        
        if len(part_inds) == 0:
            continue
        
        part_ind = mf.choice(part_inds)
        part_type = cell_part[part_ind]
        
        if part_type == 1: ## monomer with ester group
            
            process = mf.choice((mv.sci_ester, mv.sci_direct, mv.ester), p=\
                    (d1, d2, d3)/np.sum((d1, d2, d3)))
            
            if process == mv.sci_ester:
                
                cell_part[part_ind] = 0
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0].astype(int)
                cell_part[ester_ind] = -1
            
            elif process == mv.sci_direct:
                
                continue
            
            else: ## process == mv.ester
                
                cell_part[part_ind] = 0
                
                ester_ind = np.where(np.isnan(cell_part[:, 0]))[0][0].astype(int)
                cell_part[ester_ind] = -1
            
        elif part_type == 0: ## monomer W/O ester group
            
            continue
            
        elif part_type == -1: ## ester group
            
            cell_part[part_ind] = np.nan
            CH4_matrix[Z, XY, x, y, z] += 1
            
            easter_decay = mf.choice((mv.ester_CO, mv.ester_CO2), p=(d_k2, d_k3))
            
            if easter_decay == mv.ester_CO:
                CO_matrix[Z, XY, x, y, z] += 1
            else:
                CO2_matrix[Z, XY, x, y, z] += 1
            
        else: ## part_type == -2, broken ester group
            
            print('WTF???')
        

for Z, XY, x, y, z in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
    
    n_ester = len(np.where(part_matrix[Z, XY, x, y, z] == -1)[0])
    
    if n_ester != 0:
        
        ester_matrix[Z, XY, x, y, z] = n_ester

#%%
# n_H = 54 854

n_CO_0 = np.sum(CO_matrix)
n_CO2_0 = np.sum(CO2_matrix)
n_CH4_0 = np.sum(CH4_matrix)
n_ester_0 = np.sum(ester_matrix)

#%% destroy some ester groups
ester_part = 0.5

add_CO = n_ester_0 * ester_part * d_k2
n_CO = n_CO_0 + add_CO

add_CO2 = n_ester_0 * ester_part * d_k3
n_CO2 = n_CO2_0 + add_CO2

n_CH4 = n_CH4_0 + (add_CO + add_CO2)
n_ester = n_ester_0 - (add_CO + add_CO2)


print('n_H =', 54854)
print('n_ester =', n_ester)
print('n_CO =', n_CO)
print('n_CO2 =', n_CO2)
print('n_CH4 =', n_CH4)
