#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import my_constants as mc
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mc = importlib.reload(mc)

os.chdir(mv.sim_path_MAC + 'Ionization')

#%% Binding energies and occupancies
##            1s 2s 2p 3s 3p
binding_Si = [1844, 154, 104, 13.46, 8.15]
occupancy_Si = [2, 2, 6, 2, 2]

#%%
def get_gryz_ds_ddE(E, Ui, W):
    
#    inds = np.where(np.logical_and(W>Ui, W < E))
    inds = np.where(np.logical_and(W>=Ui, W <= (E+Ui)/2))
    
    dE = W[inds]
    
    diff_cs = np.pi * mc.k_el**2 * mc.e**4 / np.power(dE*mc.eV, 3) * Ui/E *\
        np.power(E / (E + Ui), 3/2) * np.power((1 - dE/E), Ui/(Ui+dE)) *\
        (dE/Ui * (1 - Ui/E) + 4/3 * np.log(2.7 + np.sqrt((E - dE)/Ui)))
    
    gryz_ds_ddE = np.zeros(len(W))
    gryz_ds_ddE[inds] = diff_cs ## m^2 / J
    
    return gryz_ds_ddE * (100)**2 * mc.eV ## cm^2 / eV


def get_total_Gryzinsky_cs(E_arr, Ui):
    
    total_Gryzinsky_cs = np.zeros(len(E_arr))
    
    dE_arr = np.logspace(0, 5, 50000)
    
    for i in range(len(E_arr)):
        
        diff_cs = get_gryz_ds_ddE(E_arr[i], Ui, dE_arr)
        total_Gryzinsky_cs[i] = np.trapz(diff_cs, x=dE_arr)
        
    return total_Gryzinsky_cs

#%% get Si diff CS
E = 4000
dE_arr = np.logspace(0, 4, 10000)

gryz_ds_ddE = np.zeros(len(dE_arr))

for i in range(len(binding_Si)):
    
    gryz_ds_ddE += occupancy_Si[i] * get_gryz_ds_ddE(E, binding_Si[i], dE_arr)

plt.loglog(dE_arr, gryz_ds_ddE * mc.n_Si, label='My')

arr = np.loadtxt('curves/4000eV.txt')

plt.loglog(arr[:, 0], arr[:, 1], label='sun hui v chai')

plt.grid()
plt.legend()
plt.show()

#%% get Si total CS
## compare wuth sun hui v chai
E_arr = np.logspace(0, 5, 1000)

total_cs_1s = get_total_Gryzinsky_cs(E_arr, binding_Si[0])*occupancy_Si[0]
total_cs_2s = get_total_Gryzinsky_cs(E_arr, binding_Si[1])*occupancy_Si[1]
total_cs_2p = get_total_Gryzinsky_cs(E_arr, binding_Si[2])*occupancy_Si[2]
total_cs_3s = get_total_Gryzinsky_cs(E_arr, binding_Si[3])*occupancy_Si[4]
total_cs_3p = get_total_Gryzinsky_cs(E_arr, binding_Si[4])*occupancy_Si[4]

plt.loglog(E_arr, total_cs_1s, label='1s')
plt.loglog(E_arr, total_cs_2s, label='2s')
plt.loglog(E_arr, total_cs_2p, label='2p')
plt.loglog(E_arr, total_cs_3s, label='3s')
plt.loglog(E_arr, total_cs_3p, label='3p')

## sun hui v chai
arr_1s = np.loadtxt('curves/1s.txt')
arr_2s = np.loadtxt('curves/2s.txt')
arr_2p = np.loadtxt('curves/2p.txt')
arr_3s = np.loadtxt('curves/3s.txt')
arr_3p = np.loadtxt('curves/3p.txt')

plt.loglog(arr_1s[:, 0], arr_1s[:, 1], label='1s sun hiu')
plt.loglog(arr_2s[:, 0], arr_2s[:, 1], label='2s sun hiu')
plt.loglog(arr_2p[:, 0], arr_2p[:, 1], label='2p sun hiu')
plt.loglog(arr_3s[:, 0], arr_3s[:, 1], label='3s sun hiu')
plt.loglog(arr_3p[:, 0], arr_3p[:, 1], label='3p sun hiu')


plt.xlim(1e+0, 1e+5)
plt.ylim(1e-25, 1e-15)

plt.title('Si subshell ionization cross sections')
plt.xlabel('E, eV')
plt.ylabel('$\sigma_{i}(E)$, cm$^2$')
plt.legend()
plt.grid()
plt.show()

#%%
plt.loglog(E_arr, total_cs_1s, label='My n=1')
plt.loglog(E_arr, total_cs_2s+total_cs_2p, label='My n=2')
plt.loglog(E_arr, total_cs_3s+total_cs_3p, label='My n=3')

arr_K = np.loadtxt('curves/EEDL_K.txt')
arr_L = np.loadtxt('curves/EEDL_L.txt')
arr_M = np.loadtxt('curves/EEDL_M.txt')

plt.loglog(arr_K[:, 0]*1e+6, arr_K[:, 1]*1e-24, label='EEDL K')
plt.loglog(arr_L[:, 0]*1e+6, arr_L[:, 1]*1e-24, label='EEDL L')
plt.loglog(arr_M[:, 0]*1e+6, arr_M[:, 1]*1e-24, label='EEDL M')

plt.xlim(1e+0, 1e+5)
plt.ylim(1e-25, 1e-15)

plt.title('Si subshell ionization cross sections')
plt.xlabel('E, eV')
plt.ylabel('$\sigma_{i}(E)$, cm$^2$')
plt.legend()
plt.grid()
plt.show()

#%%
arr_Tung = np.loadtxt('curves/Tung_Si_IMFP.txt')

total_cs = total_cs_1s+total_cs_2s+total_cs_2p+total_cs_3s+total_cs_3p

plt.loglog(E_arr, total_cs * mc.n_Si, label='My total inv IMFP')

plt.loglog(arr_Tung[:, 0], arr_Tung[:, 1] * 1e+8, label='Tung total inv IMFP')

plt.xlim(1e+0, 1e+5)
#plt.ylim(1e+4, 1e+8)

plt.title('Si subshell ionization cross sections')
plt.xlabel('E, eV')
plt.ylabel('$\sigma_{i}(E)$, cm$^2$')
plt.legend()
plt.grid()
plt.show()


