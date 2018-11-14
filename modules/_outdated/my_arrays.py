#%% Import
#import math
import numpy as np
#import sys
#from scipy import interpolate
#from random import random

#path = '/home/fedor/Yandex.Disk/Аспирантура/Моделирование/SIM_2/arrays/'
#path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/Аспирантура/Моделирование/SIM_2/arrays/'
path = '../ARRAYS/'

#%% E and theta arrays
E_arr     = np.load(path + 'E_new.npy'    )
theta_arr = np.load(path + 'theta_new.npy')
NAMES = ['H', 'C', 'O', 'Si']

#%% EEDL tables data
H_data_EEDL_full  = np.load(path + 'H_data_EEDL_full.npy' )
C_data_EEDL_full  = np.load(path + 'C_data_EEDL_full.npy' )
O_data_EEDL_full  = np.load(path + 'O_data_EEDL_full.npy' )
Si_data_EEDL_full = np.load(path + 'Si_data_EEDL_full.npy')

DATA_EEDL_FULL = [H_data_EEDL_full, C_data_EEDL_full, O_data_EEDL_full, Si_data_EEDL_full]

#%% CROSS SECTIONS
# Total
H_cs_total  = np.load(path + 'cs_all/H_cs.npy' )
C_cs_total  = np.load(path + 'cs_all/C_cs.npy' )
O_cs_total  = np.load(path + 'cs_all/O_cs.npy' )
Si_cs_total = np.load(path + 'cs_all/Si_cs.npy')
CS_TOTAL = [H_cs_total, C_cs_total, O_cs_total, Si_cs_total]

PMMA_ATOM_CS_SUM_NORM = np.zeros((len(E_arr), 4)) # H*8 + C*5 + O*2
atom_cs_sum = np.sum(H_cs_total*8, axis=1) + np.sum(C_cs_total*5, axis=1) + np.sum(O_cs_total*2, axis=1)
PMMA_ATOM_CS_SUM_NORM[:, 1] = np.sum(H_cs_total*8, axis=1)/atom_cs_sum
PMMA_ATOM_CS_SUM_NORM[:, 2] = (np.sum(H_cs_total*8, axis=1) + np.sum(C_cs_total*5, axis=1))/atom_cs_sum
PMMA_ATOM_CS_SUM_NORM[:, 3] = (np.sum(H_cs_total*8, axis=1) + np.sum(C_cs_total*5, axis=1)\
       + np.sum(O_cs_total*2, axis=1))/atom_cs_sum

# Total with sum normed to 1
H_cs_total_sum_norm  = np.load(path + 'cs_all/H_cs_sum_norm.npy' )
C_cs_total_sum_norm  = np.load(path + 'cs_all/C_cs_sum_norm.npy' )
O_cs_total_sum_norm  = np.load(path + 'cs_all/O_cs_sum_norm.npy' )
Si_cs_total_sum_norm = np.load(path + 'cs_all/Si_cs_sum_norm.npy')
CS_TOTAL_SUM_NORM = [H_cs_total_sum_norm, C_cs_total_sum_norm,\
                     O_cs_total_sum_norm, Si_cs_total_sum_norm]

# Dirrerntial cross section integral, normed to 1
H_diff_cs_int_norm  = np.load(path + 'el_diff/H_diff_int_norm.npy' )
C_diff_cs_int_norm  = np.load(path + 'el_diff/C_diff_int_norm.npy' )
O_diff_cs_int_norm  = np.load(path + 'el_diff/O_diff_int_norm.npy' )
Si_diff_cs_int_norm = np.load(path + 'el_diff/Si_diff_int_norm.npy')
DIFF_CS_INT_NORM = [H_diff_cs_int_norm, C_diff_cs_int_norm,\
                    O_diff_cs_int_norm, Si_diff_cs_int_norm]

#%% ENERGY LOSS (WITH SUBSHELLS)
H_dE_full  = np.load(path + 'dE/H_dE.npy' )
C_dE_full  = np.load(path + 'dE/C_dE.npy' )
O_dE_full  = np.load(path + 'dE/O_dE.npy' )
Si_dE_full = np.load(path + 'dE/Si_dE.npy')
DE_FULL = [H_dE_full, C_dE_full, O_dE_full, Si_dE_full]

#%% IONIZATION STUFF
# Binding energies
H_ion_Ebind  = np.load(path + 'E_bind_arr/H_E_bind_arr.npy' )
C_ion_Ebind  = np.load(path + 'E_bind_arr/C_E_bind_arr.npy' )
O_ion_Ebind  = np.load(path + 'E_bind_arr/O_E_bind_arr.npy' )
Si_ion_Ebind = np.load(path + 'E_bind_arr/Si_E_bind_arr.npy')
ION_E_BIND = [H_ion_Ebind, C_ion_Ebind, O_ion_Ebind, Si_ion_Ebind]

# 2nd electron energies
H_ion_E2nd  = np.load(path + 'E2nd/H_ion_E2nd_files.npy' )
C_ion_E2nd  = np.load(path + 'E2nd/C_ion_E2nd_files.npy' )
O_ion_E2nd  = np.load(path + 'E2nd/O_ion_E2nd_files.npy' )
Si_ion_E2nd = np.load(path + 'E2nd/Si_ion_E2nd_files.npy')
ION_E_2ND = [H_ion_E2nd, C_ion_E2nd, O_ion_E2nd, Si_ion_E2nd] # !!! UNUSED !!!

# Spectra of secondary electrons (integrals are normed to 1)
H_ion_K_spectra   = np.load(path + 'H_K_ion_spectra_int_norm.npy'  )
C_ion_K_spectra   = np.load(path + 'C_K_ion_spectra_int_norm.npy'  )
C_ion_L1_spectra  = np.load(path + 'C_L1_ion_spectra_int_norm.npy' )
C_ion_L2_spectra  = np.load(path + 'C_L2_ion_spectra_int_norm.npy' )
C_ion_L3_spectra  = np.load(path + 'C_L3_ion_spectra_int_norm.npy' )
O_ion_K_spectra   = np.load(path + 'O_K_ion_spectra_int_norm.npy'  )
O_ion_L1_spectra  = np.load(path + 'O_L1_ion_spectra_int_norm.npy' )
O_ion_L2_spectra  = np.load(path + 'O_L2_ion_spectra_int_norm.npy' )
O_ion_L3_spectra  = np.load(path + 'O_L3_ion_spectra_int_norm.npy' )
Si_ion_K_spectra  = np.load(path + 'Si_K_ion_spectra_int_norm.npy' )
Si_ion_L1_spectra = np.load(path + 'Si_L1_ion_spectra_int_norm.npy')
Si_ion_L2_spectra = np.load(path + 'Si_L2_ion_spectra_int_norm.npy')
Si_ion_L3_spectra = np.load(path + 'Si_L3_ion_spectra_int_norm.npy')
Si_ion_M1_spectra = np.load(path + 'Si_M1_ion_spectra_int_norm.npy')
Si_ion_M2_spectra = np.load(path + 'Si_M2_ion_spectra_int_norm.npy')
Si_ion_M3_spectra = np.load(path + 'Si_M3_ion_spectra_int_norm.npy')

H_ion_spectra  = [H_ion_K_spectra]
C_ion_spectra  = [C_ion_K_spectra, C_ion_L1_spectra, C_ion_L2_spectra, C_ion_L3_spectra]
O_ion_spectra  = [O_ion_K_spectra, O_ion_L1_spectra, O_ion_L2_spectra, O_ion_L3_spectra]
Si_ion_spectra = [Si_ion_K_spectra, Si_ion_L1_spectra, Si_ion_L2_spectra, Si_ion_L3_spectra,\
                  Si_ion_M1_spectra, Si_ion_M2_spectra, Si_ion_M3_spectra]
ION_SPECTRA    = [H_ion_spectra, C_ion_spectra, O_ion_spectra, Si_ion_spectra]

# Arrays to obtain 2nd electron energy
H_ion_K_E_ext   = np.load(path + 'E_H_K_ext.npy'  )
H_ion_E_ext = [H_ion_K_E_ext]

C_ion_K_E_ext   = np.load(path + 'E_C_K_ext.npy'  )
C_ion_L1_E_ext  = np.load(path + 'E_C_L1_ext.npy' )
C_ion_L2_E_ext  = np.load(path + 'E_C_L2_ext.npy' )
C_ion_L3_E_ext  = np.load(path + 'E_C_L3_ext.npy' )
C_ion_E_ext = [C_ion_K_E_ext, C_ion_L1_E_ext, C_ion_L2_E_ext, C_ion_L3_E_ext]

O_ion_K_E_ext   = np.load(path + 'E_O_K_ext.npy'  )
O_ion_L1_E_ext  = np.load(path + 'E_O_L1_ext.npy' )
O_ion_L2_E_ext  = np.load(path + 'E_O_L2_ext.npy' )
O_ion_L3_E_ext  = np.load(path + 'E_O_L3_ext.npy' )
O_ion_E_ext = [O_ion_K_E_ext, O_ion_L1_E_ext, O_ion_L2_E_ext, O_ion_L3_E_ext]

Si_ion_K_E_ext  = np.load(path + 'E_Si_K_ext.npy' )
Si_ion_L1_E_ext = np.load(path + 'E_Si_L1_ext.npy')
Si_ion_L2_E_ext = np.load(path + 'E_Si_L2_ext.npy')
Si_ion_L3_E_ext = np.load(path + 'E_Si_L3_ext.npy')
Si_ion_M1_E_ext = np.load(path + 'E_Si_M1_ext.npy')
Si_ion_M2_E_ext = np.load(path + 'E_Si_M2_ext.npy')
Si_ion_M3_E_ext = np.load(path + 'E_Si_M3_ext.npy')
Si_ion_E_ext = [Si_ion_K_E_ext, Si_ion_L1_E_ext, Si_ion_L2_E_ext, Si_ion_L3_E_ext,\
           Si_ion_M1_E_ext, Si_ion_M2_E_ext, Si_ion_M3_E_ext]

ION_E_EXT = [H_ion_E_ext, C_ion_E_ext, O_ion_E_ext, Si_ion_E_ext]
