#%% Import
import numpy as np

#path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/Study/Simulation/ARRAYS_X/'
path = '/home/fedor/Yandex.Disk/Study/Simulation/ARRAYS_X/'

atoms = ['H', 'C', 'O', 'Si']
subshells = {'H': ['K'],
             'C': ['K', 'L1', 'L2', 'L3'],
             'O': ['K', 'L1', 'L2', 'L3'],
             'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
             }
PMMA_atom_weights = np.array([8, 5, 2])

#%% E and theta arrays
E_arr = np.load(path + 'E_arr_x.npy')
theta_arr = np.load(path + 'theta_arr.npy')

#%% Cross sections
ATOMS_CS = [np.load(path + 'CS/' + atom + '_CS.npy') for atom in atoms]

#%% PMMA atoms cross sections, normed to 1
ATOM_CS_SUM = np.array([
        np.sum(ATOMS_CS[i], axis=1) for i in range(len(atoms))
        ]).transpose()

PMMA_ATOM_CS_SUM = ATOM_CS_SUM[:, :len(atoms) - 1] * PMMA_atom_weights

PMMA_ATOM_CS_SUM_NORM = np.array([
        np.sum(PMMA_ATOM_CS_SUM[:, :i], axis=1) for i in range(len(atoms))               
        ]).transpose() / np.sum(PMMA_ATOM_CS_SUM, axis=1).reshape((len(E_arr), 1))

#%% Total with sum normed to 1
ATOMS_CS_SUM_NORM = [np.array([
        np.sum(CS[:, :i], axis=1) for i in range(len(CS[0, :]) + 1)
                ]).transpose() / np.sum(CS, axis=1).reshape((len(E_arr), 1))
        for CS in ATOMS_CS]

#%%
## Dirrerntial cross section integral, normed to 1
ATOMS_DIFF_CS_INT_NORM = [np.load(path + 'DIFF_CS_int_norm/' + atom + '_DIFF_CS_int_norm.npy')
        for atom in atoms]

#%% Excitation energy loss
ATOMS_EXC_DE = [np.load(path + 'EXC_dE/' + atom + '_EXC_dE.npy') for atom in atoms]

#%% Ionization stuff
## Binding energies
ATOMS_ION_E_BIND = [np.load(path + 'ION_Ebind/' + atom + '_ION_Ebind.npy') for atom in atoms]

## 2nd electron energies
ION_E_2ND = [np.load(path + 'ION_E2nd/' + atom + '_ION_E2nd.npy') for atom in atoms]

## Spectra of secondary electrons (integrals are normed to 1)
ATOMS_ION_SPECTRA = [[np.load(path + 'ION_spectra/' + atom + '_ION_' + subshell +
        '_spectra.npy') for subshell in subshells[atom]] for atom in atoms]

## Electron energy arrays in ionization spectra
ATOMS_ION_E_SPECTRA = [[np.load(path + 'ION_Espectra/' + atom + '_ION_' + subshell +
        '_Espectra.npy') for subshell in subshells[atom]] for atom in atoms]
