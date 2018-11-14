import numpy as np
#import matplotlib.pyplot as plt
#import my_arrays as ma
#import my_functions as mf
#import my_variables as mv
#import importlib
import my_array_creator as mac
#mac = importlib.reload(mac)

#import imp
#imp.reload(mac)

#%%
arr_folder = 'IND_ARRAYS/'

el_list = [
        ['H',  'K'],
        ['C',  'K', 'L1', 'L2', 'L3'],
        ['O',  'K', 'L1', 'L2', 'L3'],
        ['Si', 'K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3']
        ]

for el in el_list:
    
    print(el)
    
    el_name = el[0]
    ss_names = el[1:]
    
    el_tot_cs = mac.get_elastic_total_cs(el_name)
    el_diff_cs = mac.get_elastic_diff_cs(el_name)
    np.save(arr_folder + el_name + '/' + el_name + '_el_tot_cs.npy', el_tot_cs)
    np.save(arr_folder + el_name + '/' + el_name + '_el_diff_cs.npy', el_diff_cs)
    
    exc_cs = mac.get_excitation_total_cs(el_name)
    exc_dE = mac.get_excitation_dE(el_name)
    np.save(arr_folder + el_name + '/' + el_name + '_exc_cs.npy', exc_cs)
    np.save(arr_folder + el_name + '/' + el_name + '_exc_dE.npy', exc_dE)
    
    
    total_cs = np.zeros(len(mac.E_arr_x))
    
    total_cs += el_tot_cs
    total_cs += exc_cs.npy
    
    for ss_name in ss_names:
        
        print(ss_name)
        
        ion_cs = mac.get_ionization_cs(el_name, ss_name)
        ion_E2nd = mac.get_ionization_E2nd(el_name, ss_name)
        ion_spectra = mac.get_ionization_2nd_spectra(el_name, ss_name, ion_E2nd)
        np.save(arr_folder + el_name + '/' + el_name + '_' + ss_name + '_ion_cs.npy', ion_cs)
        np.save(arr_folder + el_name + '/' + el_name + '_' + ss_name + '_ion_E2nd.npy', ion_E2nd)
        np.save(arr_folder + el_name + '/' + el_name + '_' + ss_name + '_ion_spectra.npy', ion_spectra)
    
        total_cs += ion_cs
    
    np.save(arr_folder + el_name + '/' + el_name + '_TOT_CS.npy', total_cs)
    