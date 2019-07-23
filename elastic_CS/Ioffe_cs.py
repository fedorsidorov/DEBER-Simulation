#%%
import numpy as np
import matplotlib.pyplot as plt
import importlib
import os

import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'Ioffe CS')

#%%
def get_Ioffe_data(Z):
    
    text_file = open('DATA/' + str(Z) + '.txt', 'r')
    lines = text_file.readlines()
    
    theta_corr = list(map(lambda s: '0'+s, lines[12].split()[2:]))
    
    theta_arr_deg_final = np.array(list(map(float, theta_corr)))
    
    E_arr_list = []
    diff_cs_list = []
    total_cs_list = []
    
    for line in lines[14:]:
        
        cs_line = line[1:]
        
        now_diff_cs, now_total_cs = cs_line.split(sep='|')
        
        now_e_diff_cs_corr = list(map(lambda s: '0'+s, now_diff_cs.split()))
        now_total_cs_corr = float('0' + now_total_cs.strip())
        
        E_arr_list.append(now_e_diff_cs_corr[0])
        diff_cs_list.append(now_e_diff_cs_corr[1:])
        total_cs_list.append(now_total_cs_corr)
    
    
    diff_cs = np.zeros((len(diff_cs_list), len(diff_cs_list[0])))
    
    for i in range(len(diff_cs_list)):
        
        diff_cs[i, :] = np.array(list(map(float, diff_cs_list[i])))
        
    
    total_cs = np.array(list(map(float, total_cs_list)))
    
    E_arr = np.array(list(map(float, E_arr_list)))
    
    E_arr_final = E_arr[::-1]
    diff_cs_final = diff_cs[::-1]
    total_cs_final = total_cs[::-1]

    return E_arr_final, theta_arr_deg_final, diff_cs_final, total_cs_final

def get_Ruth_CS(Z, E, theta):
    
    eV = 1.6e-19    
    h = 6.63e-34
    
    eps_0 = 8.85e-12
    k = 1 / (4 * np.pi * eps_0)
    
    m = 9.11e-31
    e = 1.6e-19
    
    alpha = k**2 * (m * e**4 * np.pi**2 * Z**(2/3)) / (h**2 * E * eV)
#    alpha = 0
    
    diff_cs = k**2 * Z**2 * e**4/ (4 * (E * eV)**2) / np.power(1 - np.cos(theta) + alpha, 2)
    
    return diff_cs

#%% Ar, Z = 18, E = 1000 eV
Z = 18
E = 1000

E_arr, theta_arr, diff_cs, total_cs = get_Ioffe_data(18)

diff_cs_Ruth = get_Ruth_CS(Z, E, np.deg2rad(theta_arr))

plt.semilogy(theta_arr, diff_cs[19] * 1e+16, label='Mott cross section')
plt.semilogy(theta_arr, diff_cs_Ruth * 1e+20, label='Rutherford cross section')


Ar_cs_exp = np.loadtxt('Ar_1000eV/Iga.txt')
plt.semilogy(Ar_cs_exp[:, 0], Ar_cs_exp[:, 1], 'ro', label='exp')

plt.legend()
plt.ylim(1e-3, 1e+2)

plt.xlabel(r'$\theta$ (deg)')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}$')

plt.grid()
plt.show()





