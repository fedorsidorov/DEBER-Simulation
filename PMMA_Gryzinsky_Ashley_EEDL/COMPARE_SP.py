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

#%% PMMA
E = np.logspace(0, 4.4, 1000)

SP_D = np.load('Dapor/SP_PMMA_DAPOR.npy')
SP_B = np.load('Bethe/SP_PMMA_BETHE.npy')


Aktary = np.loadtxt('Aktary.txt')
plt.loglog(Aktary[:, 0], Aktary[:, 1]*1e+8, label='Aktary')


plt.loglog(E, SP_D, label='Dapor')
plt.loglog(E, SP_B, label='Bethe')

plt.title('PMMA stopping power')
plt.xlabel('E, eV')
plt.ylabel('dE/ds, eV/cm')

plt.ylim(plt.ylim(1e+4, 1e+9))

plt.legend()
plt.grid()
plt.show()
plt.savefig('PMMA stopping power.png', dpi=300)

#%% Si
E = np.logspace(0, 4.4, 1000)
E_P = np.load('Palik/E_PALIK.npy')

SP_P = np.load('Palik/SP_SI_PALIK.npy')
SP_B = np.load('Bethe/SP_SI_BETHE.npy')

plt.loglog(E_P, SP_P, label='Palik')
plt.loglog(E, SP_B, label='Bethe')

plt.title('PMMA stopping power')
plt.xlabel('E, eV')
plt.ylabel('dE/ds, eV/cm')

plt.ylim(1e+4, 1e+9)

plt.legend()
plt.grid()
plt.show()
plt.savefig('Si stopping power.png', dpi=300)
