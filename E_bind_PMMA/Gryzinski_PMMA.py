#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import my_constants as mc
import matplotlib.pyplot as plt
import E_loss_functions as elf

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mc = importlib.reload(mc)
elf = importlib.reload(elf)

os.chdir(mv.sim_path_MAC + 'E_loss')

#%% Binding energies and occupancies
##            1s     2s     2p
binding_H  =  [13.6]
occupancy_H = [1   ]
binding_C  =  [296 , 16.59, 11.26]
occupancy_C = [2,    2,     2]
binding_O  =  [538,  28.48, 13.62]
occupancy_O = [2,    2,     4]

EE = np.logspace(0, 4.4, 1000)

#%% get total CS
CS_C_1S = elf.get_Gryzinsky_CS(EE, binding_C[0])*occupancy_C[0]
U_C_K = CS_C_1S * mc.n_PMMA*5

CS_O_1S = elf.get_Gryzinsky_CS(EE, binding_O[0])*occupancy_O[0]
U_O_K = CS_O_1S * mc.n_PMMA*2

U_PMMA_CORE = U_C_K + U_O_K

#%%
SP_C_K = elf.get_Gryzinski_SP(EE, binding_C[0], mc.n_PMMA*5, occupancy_C[0])
SP_O_K = elf.get_Gryzinski_SP(EE, binding_O[0], mc.n_PMMA*2, occupancy_O[0])

SP_PMMA_CORE = SP_C_K + SP_O_K




