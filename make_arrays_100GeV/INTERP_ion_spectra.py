import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../MODULES')
import my_functions as mf

source_folder = '../ARRAYS_100GeV/'
dest_folder = '../ARRAYS_100GeV/'

atoms = ['H', 'C', 'O', 'Si']

subshells = {'H': ['K'],
             'C': ['K', 'L1', 'L2', 'L3'],
             'O': ['K', 'L1', 'L2', 'L3'],
             'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
             }

E_arr_100GeV = np.load(source_folder + 'E_arr_100GeV.npy')

## Max value of Espectra[0] = 0.01
Ebind_arr_100GeV = np.power(10, np.arange(-2, 10, 0.01))

#%%
E_min = 100
E_max = 0

#for atom in atoms:
for atom in ['H']:
    
#    for ss in subshells[atom]:
    for ss in ['K']:
    
        atom_ION_Espectra_int_norm_raw = np.load(source_folder + 'ION_Espectra_raw/' + atom + '_ION_' +
                                    ss + '_Espectra.npy')
        atom_ION_spectra_int_norm_raw = np.load(source_folder + 'ION_spectra_raw/' + atom + '_ION_' +
                                    ss + '_spectra.npy')
        
        atom_ION_spectra_int_norm = np.zeros([len(E_arr_100GeV), len(Ebind_arr_100GeV)])
        
        for i in range(len(E_arr_100GeV)):
            
            atom_ION_spectra_int_norm[i, :] = mf.log_interp1d(atom_ION_Espectra_int_norm_raw,\
                            atom_ION_spectra_int_norm_raw[i, :])(Ebind_arr_100GeV)
            
            if i== 0:
                plt.loglog(atom_ION_Espectra_int_norm_raw,\
                           atom_ION_spectra_int_norm_raw[i, :])

#%%
for i in range(1, 1000, 50):
    plt.loglog(Ebind_arr_100GeV, atom_ION_spectra_int_norm[i, :])
