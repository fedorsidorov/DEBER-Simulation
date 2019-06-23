#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)

from scipy.optimize import curve_fit

os.chdir(mv.sim_path_MAC + 'Greeneich_curves')

#%%
def get_R(M, R0, beta, alpha):
    
    return R0 + np.power(beta**(1/alpha) / M, alpha)

#%%
data_1 = np.loadtxt('Rate/1.txt')
data_2 = np.loadtxt('Rate/2.txt')
data_3 = np.loadtxt('Rate/3.txt')
data_4 = np.loadtxt('Rate/4.txt')
data_5 = np.loadtxt('Rate/5.txt')
data_6 = np.loadtxt('Rate/6.txt')

plt.loglog(data_1[:, 0], data_1[:, 1], 'ro')
plt.loglog(data_2[:, 0], data_2[:, 1], 'ro')
plt.loglog(data_3[:, 0], data_3[:, 1], 'ro')
plt.loglog(data_4[:, 0], data_4[:, 1], 'ro')
plt.loglog(data_5[:, 0], data_5[:, 1], 'ro')
plt.loglog(data_6[:, 0], data_6[:, 1], 'ro')


M = np.logspace(2, 6, 100)

R_1 = get_R(M, 0, 9.332e+14, 3.86)
R_2 = get_R(M, 0, 1.046e+16, 3.86)

R_31 = get_R(M[:42], 0, 6.7e+9  , 2    )
R_32 = get_R(M[42:], 0, 6.645e+6, 1.188)

R_4 = get_R(M, 84, 3.14e+8 , 1.5)
R_5 = get_R(M, 241.9, 5.669e+8, 1.5)
R_6 = get_R(M, 464, 1.435e+9, 1.5)


plt.loglog(M, R_1)
plt.loglog(M, R_2)
plt.loglog(M, np.concatenate((R_31, R_32)))
plt.loglog(M, R_4)
plt.loglog(M, R_5)
plt.loglog(M, R_6)

plt.xlim(1e+3, 1e+6)
plt.ylim(1e+1, 1e+4)

plt.xlabel('M$_f$')
plt.ylabel('R, $\AA/мин$')

plt.grid()
plt.show()

plt.savefig('Greeneich_curves_cut.png', dpi=300)
