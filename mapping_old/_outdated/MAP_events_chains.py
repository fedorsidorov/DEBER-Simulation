#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
             
os.chdir(sim_path + 'mapping')

import sys
sys.path.append(sim_path + 'MODULES')

import importlib

import my_functions as mf
mf = importlib.reload(mf)

import copy

#%%
#e_data_source_dir = 'Sharma/CUBES_Sharma_2C_exc/'
e_data_source_dir = 'Sharma/CUBES_Sharma_2H_4C_exc/'
chain_source_dir = 'Sharma/CHAINS_Sharma_cube/'

## prepare e data
e_data_bank = [[[] for i in range(100)] for j in range(16)]

for k in range(16):
    for ij in range(100):
        
        e_data_bank[k][ij] = np.load(e_data_source_dir + 'cube_' + str(k) +
                   '_' + str(ij) + '.npy')

#%%
e_data_arr = copy.deepcopy(e_data_bank)

initial_L_arr = []
final_L_arr = []

out_line = np.array([100, 16, 20, 20, 20])

N_chains = 12709

n0 = 6776
n1 = n0 + 1

for n in range(N_chains):
#for n in range(n0, n1):
    
#    mf.upd_progress_bar(i, N_chains)
    if n % 1000 == 0:
        print(n, 'chains are analysed')
    
    chain_cubed_arr = np.load(chain_source_dir + 'chain_shift_' + str(n) + '_cubed.npy')
    L = len(chain_cubed_arr)
    
    initial_L_arr.append(L)
    
    scissions = []
    
    for i in range(L):
        
        line = chain_cubed_arr[i, :]
        
        if np.all(np.equal(line, out_line)):
            continue

        cube_xy = int(line[0])
        cube_z = int(line[1])
        cell_coords = line[2:]
        
        j = 0
        
        while j < len(e_data_arr[cube_z][cube_xy]):
            
            if np.all(np.equal(e_data_arr[cube_z][cube_xy][j, :], cell_coords)):

                scissions.append(i)
                e_data_arr[cube_z][cube_xy] = np.delete(e_data_arr[cube_z][cube_xy], j, axis=0)
            
            j += 1
                
    chain_edges = [0] + scissions + [L - 1]
    chain_lens = np.diff(chain_edges)

    final_L_arr += list(chain_lens)

#%%
initial_L = np.array(initial_L_arr)
final_L = np.array(final_L_arr)

for i in range(len(final_L)):
    if final_L[i] == 0:
        final_L[i] = 1

initial_log_mw = np.log10(initial_L * 100)
final_log_mw = np.log10(final_L * 100)

mf.print_histogram(initial_log_mw)
mf.print_histogram(final_log_mw)

#%%
np.save('Sharma/final_L_2H_4C_exc.npy', final_L)

#%%
initial_L = np.load('Sharma/initial_L.npy')
final_L_2C_exc = np.load('Sharma/final_L_2C_exc.npy')
final_L_4C_exc = np.load('Sharma/final_L_4C_exc.npy')
final_L_2H_4C_exc = np.load('Sharma/final_L_2H_4C_exc.npy')

for final_L in [final_L_2C_exc, final_L_4C_exc, final_L_2H_4C_exc]:
    
    for i in range(len(final_L)):
        if final_L[i] == 0:
            final_L[i] = 1

initial_log_mw = np.log10(initial_L * 100)
final_log_mw_2C_exc = np.log10(final_L_2C_exc * 100)
final_log_mw_4C_exc = np.log10(final_L_4C_exc * 100)
final_log_mw_2H_4C_exc = np.log10(final_L_2H_4C_exc * 100)
    
mf.print_histogram(initial_log_mw)
#mf.print_histogram(final_log_mw_2C_exc)
#mf.print_histogram(final_log_mw_4C_exc)
mf.print_histogram(final_log_mw_2H_4C_exc)

#%%
mat = np.loadtxt('Sharma/curves/sharma_peak_A.dat')
x_log = np.log10(mat[1:, 0])
y = mat[1:, 1]
X_LOG = np.arange(x_log[0], x_log[-1], 1e-2)
Y_LOG = mf.log_interp1d(x_log, y)(X_LOG)
#plt.plot(X_LOG, Y_LOG)
plt.semilogx(X_LOG, Y_LOG, 'ro')
plt.show()



