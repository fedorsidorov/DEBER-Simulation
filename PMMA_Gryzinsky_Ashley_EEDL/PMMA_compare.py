#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import my_constants as mc
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'PMMA_Gryzinsky_Ashley_EEDL')

#%%
## energies
E_A = np.load('Ashley/E_ASHLEY.npy')
E_D = np.load('Dapor/E_DAPOR.npy')
E_G = np.load('Gryzinsky/E_GRYZ.npy')

## total u
U_A = np.load('Ashley/U_TOTAL_ASHLEY.npy') 
U_D = np.load('Dapor/U_DAPOR.npy')
U_G = np.load('Gryzinsky/U_TOTAL_GRYZ.npy')

## core u
U_C_K_A = np.load('Ashley/U_C_K_ASHLEY.npy')
U_O_K_A = np.load('Ashley/U_O_K_ASHLEY.npy')

U_C_K_G = np.load('Gryzinsky/U_C_K_GRYZ.npy')
U_O_K_G = np.load('Gryzinsky/U_O_K_GRYZ.npy')

#%%
plt.loglog(E_A, U_A, 'o', label='Ashley total')
plt.loglog(E_A, U_C_K_A, 'o', label='Ashley C core')
plt.loglog(E_A, U_O_K_A, 'o', label='Ashley O core')

plt.loglog(E_G, U_G, '--', label='Gryzinsky total')
plt.loglog(E_G, U_C_K_G, '--', label='Gryzinsky C core')
plt.loglog(E_G, U_O_K_G, '--', label='Gryzinsky O core')

plt.loglog(E_D, U_D, label='Dapor total')

plt.title('$\mu$ for PMMA')
plt.xlabel('E, eV')
plt.ylabel('$\mu$, cm$^{-1}$')

plt.legend()
plt.grid()
plt.show()


