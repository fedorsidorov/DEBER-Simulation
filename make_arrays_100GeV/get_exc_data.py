#%%
import numpy as np
import matplotlib.pyplot as plt
import my_arrays as ma
import my_functions as mf
import my_variables as mv

E_arr_x = np.power(10, np.arange(1, 11, 0.01))

#%%
path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/Аспирантура/\
Моделирование/SIM_0/EEDL+Ioffe simulation/C_data/Excitation/'

exc_cs_raw = np.loadtxt(path + '0_83_0_0.dat', delimiter=' ', usecols=(1, 2))
exc_dE_raw = np.loadtxt(path + '0_83_11_0.dat', delimiter=' ', usecols=(1, 2))

#%%
plt.loglog(exc_dE_raw[:, 0], exc_dE_raw[:, 1])

#%%
pos_beg = 0

exc_cs = np.zeros(np.shape(E_new))
exc_dE = np.zeros(np.shape(E_new))

exc_cs[:pos_beg] = 0
exc_dE[:pos_beg] = 0

exc_cs[pos_beg:] = ms.log_interp1d(exc_cs_raw[:, 0]*1e+6,\
         exc_cs_raw[:, 1]*1e-24)(ms.E_new[pos_beg:])
exc_dE[pos_beg:] = ms.log_interp1d(exc_dE_raw[:, 0]*1e+6,\
         exc_dE_raw[:, 1]*1e+6)(ms.E_new[pos_beg:])

#%%
plt.loglog(exc_cs_raw[:, 0]*1e+6, exc_cs_raw[:, 1]*1e-24, 'ro')
plt.loglog(ms.E_new, exc_cs)
#plt.semilogx(exc_dE_raw[:, 0]*1e+6, exc_dE_raw[:, 1]*1e+6, 'ro')
#plt.semilogx(ms.E_new, exc_dE)

plt.show()

#%%
C_exc_files = np.concatenate(([exc_cs], [exc_dE]), axis=0).transpose()
