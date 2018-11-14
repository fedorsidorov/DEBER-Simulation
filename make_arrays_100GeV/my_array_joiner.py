import numpy as np
import matplotlib.pyplot as plt
#import my_arrays as ma
#import my_functions as mf
#import my_variables as mv
#import importlib
#import my_array_creator as mac
#mac = importlib.reload(mac)

arr_folder = 'IND_ARRAYS/'

#%%
E_arr = np.load(arr_folder + 'E_arr.npy')
E_arr_x = np.load(arr_folder + 'E_arr_x.npy')
theta_arr = np.load(arr_folder + 'theta_arr.npy')

H_cs_total  = np.load(arr_folder + 'H/H_TOTAL_CS.npy'  )
C_cs_total  = np.load(arr_folder + 'C/C_TOTAL_CS.npy'  )
O_cs_total  = np.load(arr_folder + 'O/O_TOTAL_CS.npy'  )
Si_cs_total = np.load(arr_folder + 'Si/Si_TOTAL_CS.npy')
CS_TOTAL = [H_cs_total, C_cs_total, O_cs_total, Si_cs_total]

PMMA_ATOM_CS_SUM_NORM = np.zeros((len(E_arr_x), 4)) # H*8 + C*5 + O*2
atom_cs_sum = np.sum(H_cs_total*8, axis=1) + np.sum(C_cs_total*5, axis=1) +\
    np.sum(O_cs_total*2, axis=1)
PMMA_ATOM_CS_SUM_NORM[:, 1] = np.sum(H_cs_total*8, axis=1)/atom_cs_sum
PMMA_ATOM_CS_SUM_NORM[:, 2] = (np.sum(H_cs_total*8, axis=1) +\
                     np.sum(C_cs_total*5, axis=1))/atom_cs_sum
PMMA_ATOM_CS_SUM_NORM[:, 3] = (np.sum(H_cs_total*8, axis=1) +\
                     np.sum(C_cs_total*5, axis=1)\
       + np.sum(O_cs_total*2, axis=1))/atom_cs_sum

#%%
H_cs_total_sum_norm  = np.load(arr_folder + 'H/H_CS_SUM_NORM.npy'  )
C_cs_total_sum_norm  = np.load(arr_folder + 'C/C_CS_SUM_NORM.npy'  )
O_cs_total_sum_norm  = np.load(arr_folder + 'O/O_CS_SUM_NORM.npy'  )
Si_cs_total_sum_norm = np.load(arr_folder + 'Si/Si_CS_SUM_NORM.npy')
CS_TOTAL_SUM_NORM = [H_cs_total_sum_norm, C_cs_total_sum_norm,\
                     O_cs_total_sum_norm, Si_cs_total_sum_norm]

# Dirrerntial cross section integral, normed to 1
H_diff_cs_int_norm  = np.load(arr_folder + 'H/H_el_diff_cs.npy'  )
C_diff_cs_int_norm  = np.load(arr_folder + 'C/C_el_diff_cs.npy'  )
O_diff_cs_int_norm  = np.load(arr_folder + 'O/O_el_diff_cs.npy'  )
Si_diff_cs_int_norm = np.load(arr_folder + 'Si/Si_el_diff_cs.npy')
DIFF_CS_INT_NORM = [H_diff_cs_int_norm, C_diff_cs_int_norm,\
                    O_diff_cs_int_norm, Si_diff_cs_int_norm]

#%% Excitation
H_exc_dE  = np.load(arr_folder + 'H/H_exc_dE.npy'  )
C_exc_dE  = np.load(arr_folder + 'C/C_exc_dE.npy'  )
O_exc_dE  = np.load(arr_folder + 'O/O_exc_dE.npy'  )
Si_exc_dE = np.load(arr_folder + 'Si/Si_exc_dE.npy')
DE_EXC = [H_exc_dE, C_exc_dE, O_exc_dE, Si_exc_dE]

#%% Ionization
H_ion_Ebind  = np.load(arr_folder + 'H/H_ION_EBIND.npy'  )
C_ion_Ebind  = np.load(arr_folder + 'C/C_ION_EBIND.npy'  )
O_ion_Ebind  = np.load(arr_folder + 'O/O_ION_EBIND.npy'  )
Si_ion_Ebind = np.load(arr_folder + 'Si/Si_ION_EBIND.npy')
ION_E_BIND = [H_ion_Ebind, C_ion_Ebind, O_ion_Ebind, Si_ion_Ebind]

# 2ndrary electron energy
H_ion_E2nd  = np.load(arr_folder + 'H/H_ION_E2ND.npy'   )
C_ion_E2nd  = np.load(arr_folder + 'C/C_ION_E2ND.npy'   )
O_ion_E2nd  = np.load(arr_folder + 'O/O_ION_E2ND.npy'   )
Si_ion_E2nd  = np.load(arr_folder + 'Si/Si_ION_E2ND.npy')
ION_E_2ND = [H_ion_E2nd, C_ion_E2nd, O_ion_E2nd, Si_ion_E2nd]

# Spectra of secondary electrons (integrals are normed to 1)
H_ion_K_spectra   = np.load(arr_folder + 'H/H_K_ion_spectra.npy'  )
C_ion_K_spectra   = np.load(arr_folder + 'C/C_K_ion_spectra.npy'  )
C_ion_L1_spectra  = np.load(arr_folder + 'C/C_L1_ion_spectra.npy' )
C_ion_L2_spectra  = np.load(arr_folder + 'C/C_L2_ion_spectra.npy' )
C_ion_L3_spectra  = np.load(arr_folder + 'C/C_L3_ion_spectra.npy' )
O_ion_K_spectra   = np.load(arr_folder + 'O/O_K_ion_spectra.npy'  )
O_ion_L1_spectra  = np.load(arr_folder + 'O/O_L1_ion_spectra.npy' )
O_ion_L2_spectra  = np.load(arr_folder + 'O/O_L2_ion_spectra.npy' )
O_ion_L3_spectra  = np.load(arr_folder + 'O/O_L3_ion_spectra.npy' )
Si_ion_K_spectra  = np.load(arr_folder + 'Si/Si_K_ion_spectra.npy' )
Si_ion_L1_spectra = np.load(arr_folder + 'Si/Si_L1_ion_spectra.npy')
Si_ion_L2_spectra = np.load(arr_folder + 'Si/Si_L2_ion_spectra.npy')
Si_ion_L3_spectra = np.load(arr_folder + 'Si/Si_L3_ion_spectra.npy')
Si_ion_M1_spectra = np.load(arr_folder + 'Si/Si_M1_ion_spectra.npy')
Si_ion_M2_spectra = np.load(arr_folder + 'Si/Si_M2_ion_spectra.npy')
Si_ion_M3_spectra = np.load(arr_folder + 'Si/Si_M3_ion_spectra.npy')

H_ion_spectra  = [H_ion_K_spectra]
C_ion_spectra  = [C_ion_K_spectra, C_ion_L1_spectra, C_ion_L2_spectra, C_ion_L3_spectra]
O_ion_spectra  = [O_ion_K_spectra, O_ion_L1_spectra, O_ion_L2_spectra, O_ion_L3_spectra]
Si_ion_spectra = [Si_ion_K_spectra, Si_ion_L1_spectra, Si_ion_L2_spectra, Si_ion_L3_spectra,\
                  Si_ion_M1_spectra, Si_ion_M2_spectra, Si_ion_M3_spectra]
ION_SPECTRA    = [H_ion_spectra, C_ion_spectra, O_ion_spectra, Si_ion_spectra]

# Arrays of extended energy to obtain 2nd electron energy
H_ion_K_E_ext   = np.load(arr_folder + 'H/H_K_Eext.npy')
H_ion_E_ext = [H_ion_K_E_ext]

C_ion_K_E_ext   = np.load(arr_folder + 'C/C_K_Eext.npy' )
C_ion_L1_E_ext  = np.load(arr_folder + 'C/C_L1_Eext.npy')
C_ion_L2_E_ext  = np.load(arr_folder + 'C/C_L2_Eext.npy')
C_ion_L3_E_ext  = np.load(arr_folder + 'C/C_L3_Eext.npy')
C_ion_E_ext = [C_ion_K_E_ext, C_ion_L1_E_ext, C_ion_L2_E_ext, C_ion_L3_E_ext]

O_ion_K_E_ext   = np.load(arr_folder + 'O/O_K_Eext.npy'  )
O_ion_L1_E_ext  = np.load(arr_folder + 'O/O_L1_Eext.npy' )
O_ion_L2_E_ext  = np.load(arr_folder + 'O/O_L2_Eext.npy' )
O_ion_L3_E_ext  = np.load(arr_folder + 'O/O_L3_Eext.npy' )
O_ion_E_ext = [O_ion_K_E_ext, O_ion_L1_E_ext, O_ion_L2_E_ext, O_ion_L3_E_ext]

Si_ion_K_E_ext  = np.load(arr_folder + 'Si/Si_K_Eext.npy' )
Si_ion_L1_E_ext = np.load(arr_folder + 'Si/Si_L1_Eext.npy')
Si_ion_L2_E_ext = np.load(arr_folder + 'Si/Si_L2_Eext.npy')
Si_ion_L3_E_ext = np.load(arr_folder + 'Si/Si_L3_Eext.npy')
Si_ion_M1_E_ext = np.load(arr_folder + 'Si/Si_M1_Eext.npy')
Si_ion_M2_E_ext = np.load(arr_folder + 'Si/Si_M2_Eext.npy')
Si_ion_M3_E_ext = np.load(arr_folder + 'Si/Si_M3_Eext.npy')
Si_ion_E_ext = [Si_ion_K_E_ext, Si_ion_L1_E_ext, Si_ion_L2_E_ext, Si_ion_L3_E_ext,\
           Si_ion_M1_E_ext, Si_ion_M2_E_ext, Si_ion_M3_E_ext]

ION_E_EXT = [H_ion_E_ext, C_ion_E_ext, O_ion_E_ext, Si_ion_E_ext]


