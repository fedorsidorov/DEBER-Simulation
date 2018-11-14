import numpy as np
import matplotlib.pyplot as plt

source_folder = '../ARRAYS_100GeV/'
dest_folder = '../ARRAYS_25keV/'
atoms = ['H', 'C', 'O', 'Si']
subshells = {'H': ['K'],
             'C': ['K', 'L1', 'L2', 'L3'],
             'O': ['K', 'L1', 'L2', 'L3'],
             'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
             }

#%% E_arr and theta_arr
E_arr_100GeV = np.load(source_folder + 'E_arr_100GeV.npy')
theta_arr_100GeV = np.load(source_folder + 'theta_arr_100GeV.npy')

end_ind = np.where(E_arr_100GeV > 25e+3)[0][1]

#%% E and theta arrays
E_arr_25keV = E_arr_100GeV[:end_ind]
theta_arr_25keV = theta_arr_100GeV[:end_ind]

np.save(dest_folder + 'E_arr_25keV.npy', E_arr_25keV)
np.save(dest_folder + 'theta_arr_25keV.npy', theta_arr_25keV)

#%% Arrays with all processes cross sections
for atom in atoms:
    
    atom_CS_100GeV = np.load(source_folder + 'CS/' + atom + '_CS.npy')
    atom_CS_25keV = atom_CS_100GeV[:end_ind, :]
    
    np.save(dest_folder + 'CS/' + atom + '_CS.npy', atom_CS_25keV)

#%% Elastic differential cross sections, integrals are normed to 1
for atom in atoms:
    
    atom_DIFF_CS_itn_norm_100GeV = np.load(source_folder + 'DIFF_CS_int_norm/' +
                                           atom + '_DIFF_CS_int_norm.npy')
    
    atom_DIFF_CS_itn_norm_25keV_pre = atom_DIFF_CS_itn_norm_100GeV[:end_ind, :]
    atom_DIFF_CS_itn_norm_25keV = atom_DIFF_CS_itn_norm_25keV_pre /\
        atom_DIFF_CS_itn_norm_25keV_pre[:, -1:]
    
    np.save(dest_folder + 'DIFF_CS_int_norm/' + atom + '_DIFF_CS_int_norm.npy',\
            atom_DIFF_CS_itn_norm_25keV)

#%% Excitation dE
for atom in atoms:
    
    atom_EXC_dE_100GeV = np.load(source_folder + 'EXC_dE/' + atom + '_EXC_dE.npy')
    atom_EXC_dE_25keV = atom_EXC_dE_100GeV[:end_ind]
    
    np.save(dest_folder + 'EXC_dE/' + atom + '_EXC_dE.npy', atom_EXC_dE_25keV)

#%% Ionization Ebind
for atom in atoms:
    
    atom_ION_Ebind_100GeV = np.load(source_folder + 'ION_Ebind/' + atom + '_ION_Ebind.npy')
    atom_ION_Ebind_25keV = atom_ION_Ebind_100GeV
    
    np.save(dest_folder + 'ION_Ebind/' + atom + '_ION_Ebind.npy', atom_ION_Ebind_25keV)

#%% Ionization E2nd
for atom in atoms:
    
    atom_ION_E2nd_100GeV = np.load(source_folder + 'ION_E2nd/' + atom + '_ION_E2nd.npy')
    atom_ION_E2nd_25keV = atom_ION_E2nd_100GeV[:end_ind, :]
    
    np.save(dest_folder + 'ION_E2nd/' + atom + '_ION_E2nd.npy', atom_ION_E2nd_25keV)

#%% Ionization spectra
for atom in atoms:
        for ss in subshells[atom]:
            
            atom_ION_Espectra_100GeV = np.load(source_folder + 'ION_Espectra/' +
                                               atom + '_ION_' + ss + '_Espectra.npy')
            atom_ION_spectra_int_norm_100GeV = np.load(source_folder + 'ION_spectra/' +
                                                       atom + '_ION_' + ss + '_spectra.npy')
            
            atom_ION_ss_Espectra_25keV = atom_ION_Espectra_100GeV
            atom_ION_ss_spectra_int_norm_25keV = atom_ION_spectra_int_norm_100GeV[:end_ind, :]
            
            np.save(dest_folder + 'ION_Espectra/' + atom + '_ION_' + ss + '_Espectra.npy',
                    atom_ION_ss_Espectra_25keV)
            np.save(dest_folder + 'ION_spectra/' + atom + '_ION_' + ss + '_spectra.npy',
                    atom_ION_ss_spectra_int_norm_25keV)




