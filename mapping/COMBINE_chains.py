#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
from random import random

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
             
os.chdir(sim_path + 'mapping')

import sys
sys.path.append(sim_path + 'MODULES')

import importlib

import my_functions as mf
mf = importlib.reload(mf)

#%%
source_dir = '../make_chains/CHAINS_Sharma/'
dest_dir = 'Sharma/CHAINS_Sharma/'

N_chains = 10000

chain_bank = [[] for i in range(N_chains)]
L_arr = np.zeros(N_chains)

## load chains into bank
for i in range(N_chains):
    
    mf.upd_progress_bar(i, N_chains)
    
    chain_bank[i] = np.load(source_dir + 'chain_' + str(i % len(chain_bank)) + '.npy')
    L_arr[i] = len(chain_bank[i])
    
## check log_mw
log_mw = np.log10(L_arr * 100)

#%% prepare histograms
lx, ly, lz = 100, 100, 160
space = 100

x_beg, y_beg, z_beg = 0, 0, 0
x_end, y_end, z_end = x_beg + lx, y_beg + ly, z_beg + lz

step_prec = 10
step_mono = 0.28

x_bins_total = np.array([x_beg, x_end])
y_bins_total = np.array([y_beg, y_end])
z_bins_total = np.array([z_beg, z_end])

x_bins_prec = np.arange(x_beg, x_end + 1, step_prec)
y_bins_prec = np.arange(y_beg, y_end + 1, step_prec)
z_bins_prec = np.arange(z_beg, z_end + 1, step_prec)

#x_bins_mono = np.arange(x_beg, x_end/2, step_mono)
#y_bins_mono = np.arange(y_beg, y_end/2, step_mono)
#z_bins_mono = np.arange(z_beg, z_end/2, step_mono)

bins_total = [x_bins_total, y_bins_total, z_bins_total]
bins_prec = [x_bins_prec, y_bins_prec, z_bins_prec]
#bins_mono = [x_bins_mono, y_bins_mono, z_bins_mono]

hist_total = np.zeros((len(x_bins_total) - 1, len(y_bins_total) - 1, len(z_bins_total) - 1))
hist_prec = np.zeros((len(x_bins_prec) - 1, len(y_bins_prec) - 1, len(z_bins_prec) - 1))
#hist_mono = np.zeros((len(x_bins_mono) - 1, len(y_bins_mono) - 1, len(z_bins_mono) - 1))

#%% create chain_list and check density
chain_list = []

V = lx*ly*lz * (1e-7)**3 ## cm^3
m_mon = 1.66e-22 ## g
rho = 1.19 ## g / cm^3

i = 0

while True:
    
    if i % 100 == 0:
        print(i, 'chains are added')
    
    if np.sum(hist_total) * m_mon / V >= rho:
        print('Needed density is achieved')
        break

    else:
        
        ii = i % len(chain_bank)
        new_chain = chain_bank[ii]
    
        x_shift = mf.uniform(-space, lx + space)
        y_shift = mf.uniform(-space, ly + space)
        
        new_chain_shift = new_chain + np.array((x_shift, y_shift, 0))
        
        if new_chain_shift.max(axis=0)[0] < x_beg or new_chain_shift.max(axis=0)[1] < y_beg\
        or new_chain_shift.min(axis=0)[0] > x_end or new_chain_shift.min(axis=0)[1] > y_end:
            continue
        
        chain_list.append(new_chain_shift)
        
        hist_total += np.histogramdd(new_chain_shift, bins=bins_total)[0]
        hist_prec += np.histogramdd(new_chain_shift, bins=bins_prec)[0]
#        hist_mono += np.histogramdd(new_chain_shift, bins=bins_mono)[0]

    i += 1

#%% check density
density_total = hist_total[0][0][0] * m_mon / V
density_precise = hist_prec * m_mon / V * (len(x_bins_prec) - 1) * (len(y_bins_prec) - 1) *\
    (len(z_bins_prec) - 1)


#%% save chains to files
i = 0

for chain in chain_list:
    
    mf.upd_progress_bar(i, len(chain_list))
    np.save(dest_dir + 'chain_shift_' + str(i) + '.npy', chain)
    i += 1

#%% check n_mon in 0.28 nm cube
#array_mono = np.reshape(hist_mono, (np.prod(np.shape(hist_mono)),))
#array_mono_hist = np.histogram(array_mono, bins=10)[0]
#mf.print_histogram(array_mono, user_bins=10, is_normed=False)

#%% cut chains to cube shape
#chain_cut_list = []
#
#for chain in chain_list:
#    
#    statements = [chain[:, 0] >= x_min, chain[:, 0] <= x_max,
#                  chain[:, 1] >= x_min, chain[:, 1] <= x_max]
#    inds = np.where(np.logical_and.reduce(statements))[0]
#    
#    beg = 0
#    end = -1
#    
#    for i in range(len(inds) - 1):
#        if inds[i+1] > inds[i] + 1 or i == len(inds) - 2:
#            end = i + 1
#            chain_cut_list.append(chain[inds[beg:end], :])
#            beg = i + 1

#%% get nice 3D picture
#l_xy = len(np.arange(100, 151, 1))
#l_z = len(np.arange(0, 81, 1))
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
#for chain in chain_list[0:-1:5]:
#    ax.plot(chain[:, 0], chain[:, 1], chain[:, 2])
#
#ax.plot(np.arange(100, 151, 1), np.ones(l_xy)*100, np.ones(l_xy)*0, 'k')
#ax.plot(np.arange(100, 151, 1), np.ones(l_xy)*150, np.ones(l_xy)*0, 'k')
#ax.plot(np.arange(100, 151, 1), np.ones(l_xy)*100, np.ones(l_xy)*80, 'k')
#ax.plot(np.arange(100, 151, 1), np.ones(l_xy)*150, np.ones(l_xy)*80, 'k')
#
#ax.plot(np.ones(l_xy)*100, np.arange(100, 151, 1), np.ones(l_xy)*0, 'k')
#ax.plot(np.ones(l_xy)*150, np.arange(100, 151, 1), np.ones(l_xy)*0, 'k')
#ax.plot(np.ones(l_xy)*100, np.arange(100, 151, 1), np.ones(l_xy)*80, 'k')
#ax.plot(np.ones(l_xy)*150, np.arange(100, 151, 1), np.ones(l_xy)*80, 'k')
#
#ax.plot(np.ones(l_z)*100, np.ones(l_z)*100, np.arange(0, 81, 1), 'k')
#ax.plot(np.ones(l_z)*150, np.ones(l_z)*100, np.arange(0, 81, 1), 'k')
#ax.plot(np.ones(l_z)*100, np.ones(l_z)*150, np.arange(0, 81, 1), 'k')
#ax.plot(np.ones(l_z)*150, np.ones(l_z)*150, np.arange(0, 81, 1), 'k')
#
#plt.xlim(x_beg, x_end)
#plt.ylim(y_beg, y_end)
#plt.title('Polymer chain simulation')
#ax.set_xlabel('x, nm')
#ax.set_ylabel('y, nm')
#ax.set_zlabel('z, nm')
#plt.show()
