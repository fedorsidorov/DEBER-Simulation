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
part_matrix_master = np.load('Sharma/MATRIX_particles.npy')
chain_inv_matrix = np.load('Sharma/MATRIX_chains_inv.npy')

#%%
part_matrix = copy.deepcopy(part_matrix_master) ## particles to interact with
CO2_CH4_matrix = np.zeros(np.shape(e_matrix))

#%%
s0, s1, s2, s3, s4 = np.shape(e_matrix)

chain_lens = np.load('Sharma/initial_L.npy')

bond_scissions = [[] for i in range(mv.N_chains)]

## probabilities
p1 = 0.3
p2 = 0.5
p3 = 0.7
probs = [p1 - 0, p2 - p1, p3 - p2, 1 - p3]

for Z, XY, z, y, x in product(range(s0), range(s1), range(s2), range(s3), range(s4)):
    
    if XY == 0:
        print(Z)
    
    n_events = e_matrix[Z, XY, x, y, z].astype(int)
    
    for i in range(n_events):
        
        ## process identificatior
        proc_ind = mf.choice(mv.proc_indexes, p=probs)
        
        ## only ... of 5 C atoms are of interest
        if proc_ind == mv.nothing:
            continue
        
        ## current cell of working array 
        cell_part = part_matrix[Z, XY, z, y, x]
        
        ## indexes of all particles we have - monomers and ester groups
        part_inds = np.where(np.logical_not(np.isnan(cell_part[:, 0])))[0]
        
        if len(part_inds) == 0:
            continue
        
        now_part_line = cell_part[mf.choice(part_inds), :]
        
        if np.all(now_part_line == -1): ## now particle is ester group
        
            ## remove ester group from the matrix
            cell_part[now_part_ind] = np.nan
            ## add one CO2_CH4 entry
            CO2_CH4_matrix[[Z, XY, x, y, z]] += 1
        
        else: ## now particle is monomer
            
            n_chain, n_mon, mon_type = now_part_line
            
            
            if mon_type == 0: ## bonded monomer
                
                new_mon_type = mf.choice((-1, 1))
                
                n_next_mon = n_mon + new_mon_type
                
                Z_new, XY_new, x_new, y_new, z_new, n_mon_new, pos_new =\
                    chain_inv_matrix[n_chain, ]
                
                
                if proc_ind == mv.sci_ester: ## scission with ester group
                
                elif proc_ind == mv.sci_direct: ## scission without ester group
                
                elif proc_ind == mv.ester: ## scission without ester group
            
            
            elif np.abs(mon_type) == 1: ## half bonded monomer
                        
            else: ## free monomer            
                        
                        
                        
                        
                        
                        
                        