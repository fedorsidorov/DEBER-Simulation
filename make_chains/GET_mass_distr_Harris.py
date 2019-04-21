#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os

sim_path = '/Users/fedor/Documents/DEBER-Simulation/'             
os.chdir(sim_path + 'make_chains')

import importlib

import my_functions as mf
mf = importlib.reload(mf)

#%%
mat = np.loadtxt('harris1973_A.dat')

## y is cumulative!
x = mat[:, 0]
y_c = mat[:, 1]

y_c[np.where(y_c < 0)] = 0
y_c[np.where(y_c > 1)] = 1

y_c[0] = 0
y_c[-3:] = 1

#plt.semilogx(x, y)

X_LOG = np.arange(np.log10(x[0]), np.log10(x[-1]), 1e-2)
Y_LOG_C = mf.log_interp1d(np.log10(x), y_c)(X_LOG)

X_LOG = X_LOG[np.where(np.logical_not(np.isnan(Y_LOG_C)))]
Y_LOG_C = Y_LOG_C[np.where(np.logical_not(np.isnan(Y_LOG_C)))]

plt.plot(X_LOG, Y_LOG_C, 'ro')
plt.grid()
plt.show()

#%%
def get_log_mw():
    r = mf.random()
    for i in range(len(Y_LOG_C) - 1):
        if r < Y_LOG_C[i + 1]:
            return X_LOG[i]

N_chains = 10000

log_mw_arr = np.zeros(N_chains)
L_arr = np.zeros(N_chains)

for i in range(N_chains):
    log_mw = get_log_mw()
    log_mw_arr[i] = log_mw
    L_arr[i] = int(10**log_mw / 100)

#%%
L_arr_cut = L_arr[np.where(L_arr < 5e+4)]

plt.hist(np.log10(L_arr), bins=20)
plt.hist(np.log10(L_arr_cut), bins=20)

