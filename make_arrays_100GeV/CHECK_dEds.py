import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('../MODULES')
import my_variables as mv
import my_functions as mf
import importlib
#mv = importlib.reload(mv)


#%%
path = 'ARRAYS_X'

E_arr_x = np.load(path + os.sep + 'E_arr_x.npy')
theta_arr = np.load(path + os.sep + 'theta_arr.npy')

ss_dict = {'H': ['K'],
           'C': ['K', 'L1', 'L2', 'L3'],
           'O': ['K', 'L1', 'L2', 'L3'],
           'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
        }

#%%
elem = 'Si'
full_path = path + os.sep + elem + os.sep + elem

TOTAL_CS = np.load(full_path + '_TOTAL_CS.npy')
exc_dE = np.load(full_path + '_exc_dE.npy')

l_exc_inv = mv.CONC_at[elem]*TOTAL_CS[:, 1]
dEds_exc = exc_dE*l_exc_inv

ION_EBIND = np.load(full_path + '_ION_EBIND.npy')
ION_E2ND = np.load(full_path + '_ION_E2ND.npy')

l_ion_inv = mv.CONC_at[elem]*TOTAL_CS[:, 2:]
dEds_ion = (ION_E2ND + ION_EBIND)*l_ion_inv

EEDL_dEds = np.loadtxt('.' + os.sep + 'curves' + os.sep + elem + '_dEds.txt')

plt.loglog(E_arr_x, dEds_exc, label='exc')
plt.loglog(E_arr_x, dEds_ion, label='ion')
plt.loglog(E_arr_x, (dEds_exc + np.sum(dEds_ion, axis=1)), 'ro', label='TOTAL')
plt.loglog(EEDL_dEds[:, 0], EEDL_dEds[:, 1], 'b--', label='EEDL')
plt.title(elem)
plt.xlabel('E, eV')
plt.ylabel('dE/ds, eV/cm')
plt.legend()
plt.show()





