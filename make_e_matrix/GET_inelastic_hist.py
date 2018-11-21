#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'make_e_matrix')

#%%
n_files = 1922

n_bins = 201
x_bins = np.linspace(-100, 100, n_bins)
x_arr = (x_bins[:-1] + x_bins[1:]) / 2
x_hist = np.zeros(len(x_arr))

for i in range(n_files):
    
    mf.upd_progress_bar(i, n_files)
    
    now_DATA_Pn = np.load(mv.sim_path_MAC + 'e_DATA/DATA_Pn_20keV_122nm/DATA_Pn_' +\
                       str(i) +'.npy')
    x_hist += np.histogram(now_DATA_Pn[:, 5], bins=x_bins)[0]

#%%
plt.semilogy(x_arr, x_hist / x_hist.max(), label = 'inelastic events')
plt.plot(x_arr, np.ones(np.shape(x_arr)) * 0.01, label = '1% level')
plt.xlim((-50, 50))
plt.grid()
plt.legend()
plt.show()

#%% 20 nm is ok, but 50 nm is better