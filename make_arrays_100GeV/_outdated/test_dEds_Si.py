#%% Import
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../MODULES')
#import my_arrays as ma
#import my_functions as mf
import my_array_creator as mac

#%% Second approach
dsdE_arr = np.zeros((len(mac.E_arr_x), 2))

for fn in range(500):
    
    print(fn)
#    DATA = np.load('./DATA/Si_exc/DATA_Si_' + str(fn) + '.npy')
    DATA = np.load('./DATA_PMMA_Si/DATA_PMMA_Si_' + str(fn) + '.npy')
    
    for i in range(len(DATA) - 1):
        
        if DATA[i+1, 0] == DATA[i, 0]:
            
            E = DATA[i, 4]
            dE = DATA[i+1, 8]
            
#            if dE >= 2e+4:
#                print('ERROR!')
#                continue
#            else:
            ds = np.linalg.norm((DATA[i+1, 5:8] - DATA[i, 5:8]))
#                dsdE = ds / dE
                
            E_pos = np.argmin(np.abs(mac.E_arr_x - E))
                
            dsdE_arr[E_pos, 0] += ds
            dsdE_arr[E_pos, 1] += dE

dsdE_avg = np.zeros(np.size(mac.E_arr_x))

for i in range(len(dsdE_avg)):
 
    if dsdE_arr[i, 1] == 0:
        dsdE_avg[i] = np.inf
    else:
        dsdE_avg[i] = dsdE_arr[i, 0] / dsdE_arr[i, 1]


#%%
#Si_EEDL = np.load('./arrays/Si_data_EEDL_full.npy')
dsdE_avg = np.load('./dsdE_avg.npy')
#plt.loglog(ma.E_arr, Si_EEDL[:, 7])

mat = np.loadtxt('dEds_paper.csv')
x = mat[:, 0]
y = mat[:, 1]
#y[1] = y[0]
#x = 10**x_log
#y = 10**y_log
#X = np.arange(x[0], x[-1], 1e+4)
#Y = mf.log_interp1d(x, y)(X)

mat_nist = np.loadtxt('dEds_NIST.csv')
x_nist = mat_nist[:, 0]
y_nist = mat_nist[:, 1]

plt.loglog(x, y*1e+8, linewidth=5, label='paper')
plt.loglog(x_nist*1e+6, y_nist*1e+6*1.18, linewidth=5, label='NIST')
plt.loglog(mac.E_arr_x[:500], 1/dsdE_avg[:500], 'ro', label='we')
plt.xlabel('E, eV')
plt.ylabel('dE/ds, eV/cm')
plt.grid()
plt.title('Energy loss')
plt.legend()
plt.show()




