import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import importlib

sys.path.append('../MODULES')
import my_functions as mf
import my_arrays_x_smart as ma

#mv = importlib.reload(mv)

#%%
ss_dict = {'H': ['K'],
           'C': ['K', 'L1', 'L2', 'L3'],
           'O': ['K', 'L1', 'L2', 'L3'],
           'Si': ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
        }

atoms = list(ss_dict.keys())

SPECTRA = ma.ATOMS_ION_SPECTRA
E_SPECTRA = ma.ATOMS_ION_E_SPECTRA

#%%
atom_id = 0
ss_id = 0

spectra = SPECTRA[atom_id][ss_id]
e_spectra = E_SPECTRA[atom_id][ss_id]

plt.loglog(e_spectra, spectra[20:50, :].transpose())











