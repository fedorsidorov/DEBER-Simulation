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

os.chdir(mv.sim_path_MAC + 'development')

#%%
dE_matrix = np.load('MATRIX_dE_Aktary_100uC_lines.npy')

e_loss_matrix_2 = np.sum(dE_matrix, axis=1) / 100

#%%
Mn = 950e+3
Gs = 1
rho = 1.19e-21
Na = 6.02e+23

Mf_matrix = Mn / (1 + Gs * e_loss_matrix_2 * Mn / (rho * Na))

R0 = 0.0 # nm/s
alpha = 3.86
beta = 9.332e+14

R_matrix = R0 + beta / (np.power(Mf_matrix, alpha))

T_matrix = 1 / R_matrix

profile = copy.deepcopy(T_matrix)


def dissolve(profile, T):
    
    profile -= T
    
    profile[np.where(profile < 0)] = 0
    profile[np.where(profile > 0)] = 1


T_new = dissolve(profile, 1e-7)


fig, ax = plt.subplots(1,1)

ax.imshow(profile.transpose()*(-1), cmap='Greys')

ax.set_xlabel('x, nm')
ax.set_ylabel('z, nm')

ax.set_xticklabels(['', 0, 20, 40, 60, 80, 100])
ax.set_yticklabels(['', 0, 20, 40, 60, 80, 100])
#fig.colorbar()

ax.set_title('Direct Monte-Carlo simulation')