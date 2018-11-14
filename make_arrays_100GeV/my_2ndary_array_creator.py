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

#%%
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


    # Make total_cs array - OK
    total_cs = np.zeros((len(E_arr_x), len(el) + 1))
    
    total_cs[:, 0] = np.load(arr_folder + el_name + '/' + el_name + '_el_tot_cs.npy')
    total_cs[:, 1] = np.load(arr_folder + el_name + '/' + el_name + '_exc_cs.npy')
    
    i = 2
    
    for ss_name in ss_names:
        
        total_cs[:, i] = np.load(arr_folder + el_name + '/' + el_name + '_' +\
                ss_name + '_ion_cs.npy')
        i += 1
    
#    np.save(arr_folder + el_name + '/' + el_name + '_TOTAL_CS.npy', total_cs)
    
    
    # Make total_cs_sum_norm array - OK
#    total_cs_sum_norm = np.zeros((len(E_arr_x), len(el) + 2))
#    
#    cs_sum = np.sum(total_cs, axis=1)
#    
#    for j in range(len(total_cs[0, :])):
#            
#        total_cs_sum_norm[:, j] = np.sum(total_cs[:, :j], axis=1) / cs_sum
#    
#    total_cs_sum_norm[:, -1] = 1
#    
#    np.save(arr_folder + el_name + '/' + el_name + '_CS_SUM_NORM.npy', total_cs_sum_norm)
    
    # Make ion_E2nd arrays - OK
    ion_E2nd = np.zeros((len(E_arr_x), len(ss_names)))
    
    j = 0
    
    for ss_name in ss_names:
        
        ion_E2nd[:, j] = np.load(arr_folder + el_name + '/' + el_name + '_' +\
                ss_name + '_ion_E2nd.npy')
        j += 1
    
    np.save(arr_folder + el_name + '/' + el_name + '_ION_E2ND.npy', ion_E2nd)
    
    
    