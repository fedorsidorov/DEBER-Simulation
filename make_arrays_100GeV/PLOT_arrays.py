import numpy as np
import matplotlib.pyplot as plt
import os

#%% Load E and theta arrays
path = 'IND_ARRAYS'

E_arr_x = np.load(path + os.sep + 'E_arr_x.npy')
theta_arr = np.load(path + os.sep + 'theta_arr.npy')

#%% Load individual data
elem = 'H'

ss_dict = {'H': ['K'],
           'C': ['K', 'L1', 'L2', 'L3'],
           'O': ['K', 'L1', 'L2', 'L3'],
           'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
        }

full_path = path + os.sep + elem + os.sep + elem

CS_SUM_NORM = np.load(full_path + '_CS_SUM_NORM.npy')
TOTAL_CS = np.load(full_path + '_TOTAL_CS.npy')
ION_E_BIND = np.load(full_path + '_ION_EBIND.npy')
ION_E2ND = np.load(full_path + '_ION_E2ND.npy')

el_diff_cs = np.load(full_path + '_el_diff_cs.npy')
el_tot_cs = np.load(full_path + '_el_tot_cs.npy')
exc_cs = np.load(full_path + '_exc_cs.npy')
exc_dE = np.load(full_path + '_exc_dE.npy')

for ss_name in ss_dict[elem]:
    exec(ss_name + '_Eext     = np.load(full_path + \'_' + ss_name + '_Eext.npy\')')
    exec(ss_name + '_ion_cs   = np.load(full_path + \'_' + ss_name + '_ion_cs.npy\')')
    exec(ss_name + '_ion_E2nd = np.load(full_path + \'_' + ss_name + '_ion_E2nd.npy\')')
    exec(ss_name + '_ion_spectra = np.load(full_path + \'_' + ss_name + '_ion_spectra.npy\')')

#%% Plot data
#plt.loglog(E_arr_x, TOTAL_CS)
#plt.loglog(theta_arr, el_diff_cs[250, :])








