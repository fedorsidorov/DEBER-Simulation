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

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mi = importlib.reload(mi)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'mapping')

import my_mapping as mm
mm = importlib.reload(mm)

#%%
e_matrix  = np.load(mv.sim_path_MAC + 'MATRIXES/Aktary/MATRIX_Aktary_100uC_C_ion.npy')
dE_matrix = np.load(mv.sim_path_MAC + 'MATRIXES/Aktary/MATRIX_Aktary_100uC_dE.npy')

resist_matrix = np.load(mv.sim_path_MAC + 'MATRIXES/Aktary/MATRIX_resist_Aktary.npy')
chain_table   = np.load(mv.sim_path_MAC + 'MATRIXES/Aktary/TABLE_chains_Aktary.npy')

scission_matrix = np.zeros(np.shape(e_matrix))

N_chains_total = 9331
N_mon_chain_max = 6590

sci_per_mol_matrix = np.zeros(N_chains_total)

resist_shape = np.shape(resist_matrix)[:3]

#%%
e_loss_matrix_2 = np.sum(dE_matrix / 2**3, axis=1) / 50

Mf_matrix = 950e+3 / (1 + 0.9 * e_loss_matrix_2 * 950e+3 / (1.19e-21 * 6.02e+23))

R0 = 0.0 # nm/s
alpha = 3.86
beta = 9.332e+14

R_matrix = R0 + beta / (np.power(Mf_matrix, alpha))

T_matrix = 2 / R_matrix

def dissolve(T_matrix, T):
    
#    T_half = (T_matrix[:25, :] + T_matrix[25:, :][::-1, :]) / 2
    T_half = T_matrix
    
    T_half -= T
    
    T_half[np.where(T_half < 0)] = 0
    T_half[np.where(T_half > 0)] = 1
    
#    T_new = np.vstack((T_half, T_half[::-1, :]))
    T_new = T_half
    
    return T_new
    
T_new = dissolve(T_matrix, 1e-3)

plt.imshow(T_new.transpose())
