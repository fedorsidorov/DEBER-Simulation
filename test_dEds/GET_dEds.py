#%% Import
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append('../MODULES')
import my_arrays_x as ma
import my_variables as mv

#%% Second approach
SP_arr = np.zeros((len(ma.E_arr), 2))
SP_no_el_arr = np.zeros((len(ma.E_arr), 2))
IMFP_arr = np.zeros((len(ma.E_arr), 2))

filenames = os.listdir('DATA')

for fn in range(43):
    
    print(fn)
    
    DATA = np.load('DATA/DATA_PMMA_' + str(fn) + '.npy')
    
    for i in range(len(DATA) - 1):
        
        need_new_E_pos = True
        E_pos = -1
        
        if DATA[i+1, 0] == DATA[i, 0]:
            
            if DATA[i, 3] == 0 and E_pos != -1: ## elastic scattering not in track start
                need_new_E_pos = False
            else:
                need_new_E_pos = True
            
            E = DATA[i, 4]
            dE = DATA[i, 8]
            
            if dE > 25e+3:
                print('ERROR!')
                continue
            else:
                ds = np.linalg.norm((DATA[i+1, 5:8] - DATA[i, 5:8]))
            
            if need_new_E_pos:
                E_pos = np.argmin(np.abs(ma.E_arr - E))
                
            SP_arr[E_pos, 0] += ds
            SP_arr[E_pos, 1] += dE

            if dE > 0:
                IMFP_arr[E_pos, 0] += ds
                IMFP_arr[E_pos, 1] += 1
        
        else:
            need_new_E_pos = True
            E_pos = -1

SP_avg = np.zeros(np.size(ma.E_arr))
IMFP_avg = np.zeros(np.size(ma.E_arr))

for i in range(len(SP_avg)):
 
    if SP_arr[i, 1] == 0:
        SP_avg[i] = np.inf
    else:
        SP_avg[i] = SP_arr[i, 0] / SP_arr[i, 1]
    
    if IMFP_arr[i, 1] == 0:
        IMFP_avg[i] = 0
    else:
        IMFP_avg[i] = IMFP_arr[i, 0] / IMFP_arr[i, 1]

#%% SP
SP_Tahir = np.loadtxt('curves/SP.txt')
SP_NIST = np.loadtxt('curves/SP_NIST.txt')
plt.loglog(ma.E_arr, 1/SP_avg, 'ro', label='me')
plt.loglog(SP_Tahir[:, 0], SP_Tahir[:, 1]*1e+8, '-b', label='Tahir')
plt.loglog(SP_NIST[:, 0]*1e+6, SP_NIST[:, 1]*1e+6*mv.rho_PMMA, '-b', label='NIST')
plt.xlim(1e+1, 3e+4)
plt.legend()
plt.show()

#%% IMFP
IMFP_Tahir = np.loadtxt('curves/IMFP.txt')
plt.loglog(ma.E_arr, IMFP_avg, 'ro', label='me')
plt.loglog(IMFP_Tahir[:, 0], IMFP_Tahir[:, 1]*1e-8, '-b', label='Tahir')
plt.xlim(1e+1, 3e+4)
plt.legend()
plt.show()
