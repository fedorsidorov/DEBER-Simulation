#%% Import
import numpy as np
import os
import importlib
import my_arrays_25keV as ma
import my_functions as mf
import my_variables as mv
import my_constants as mc
import matplotlib.pyplot as plt

ma = importlib.reload(ma)
mf = importlib.reload(mf)
mv = importlib.reload(mv)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'PMMA_Gryzinsky_Ashley_EEDL')

#%%
E = ma.E_arr
ATOMS_CS = ma.ATOMS_CS

U_H = np.sum(ATOMS_CS[0][:, 1:], axis=1) * mc.n_PMMA*8
U_C = np.sum(ATOMS_CS[1][:, 1:], axis=1) * mc.n_PMMA*5
U_O = np.sum(ATOMS_CS[2][:, 1:], axis=1) * mc.n_PMMA*2

U_TOTAL = U_H + U_C + U_O
