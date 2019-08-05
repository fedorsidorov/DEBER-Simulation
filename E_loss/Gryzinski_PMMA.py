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

os.chdir(mv.sim_path_MAC + 'E_loss')

#%% Binding energies and occupancies
##            1s 2s 2p 3s 3p
binding_H  =  [13.6]
occupancy_H = [1]
binding_C  =  [296, 16.59, 11.26]
occupancy_C = [2, 2, 2]
binding_O  =  [538, 28.48, 13.62]
occupancy_O = [2, 2, 4]

#%%
def get_gryz_ds_ddE(E, Ui, W):

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


#%% get total CS
E_arr = np.logspace(0, 5, 1000)

## H
CS_H_1S = get_total_Gryzinsky_cs(E_arr, binding_H[0])*occupancy_H[0]

CS_H_TOTAL = CS_H_1S
U_H_TOTAL = CS_H_TOTAL * mc.n_PMMA*8

## C
CS_C_1S = get_total_Gryzinsky_cs(E_arr, binding_C[0])*occupancy_C[0]
CS_C_2S = get_total_Gryzinsky_cs(E_arr, binding_C[1])*occupancy_C[1]
CS_C_2P = get_total_Gryzinsky_cs(E_arr, binding_C[2])*occupancy_C[2]

U_C_K = CS_C_1S * mc.n_PMMA*5
U_C_TOTAL = (CS_C_1S + CS_C_2S + CS_C_2P) * mc.n_PMMA*5

## O
CS_O_1S = get_total_Gryzinsky_cs(E_arr, binding_O[0])*occupancy_O[0]
CS_O_2S = get_total_Gryzinsky_cs(E_arr, binding_O[1])*occupancy_O[1]
CS_O_2P = get_total_Gryzinsky_cs(E_arr, binding_O[2])*occupancy_O[2]

U_O_K = CS_O_1S * mc.n_PMMA*2
U_O_TOTAL = (CS_O_1S + CS_O_2S + CS_O_2P) * mc.n_PMMA*2

#%%
U_TOTAL = U_H_TOTAL + U_C_TOTAL + U_O_TOTAL

#%%
plt.loglog(E_arr, U_H_TOTAL, label='H total')
plt.loglog(E_arr, U_C_TOTAL, label='C total')
plt.loglog(E_arr, U_O_TOTAL, label='O total')

plt.loglog(E_arr, U_H_TOTAL+U_C_TOTAL+U_O_TOTAL, label='MMA')

plt.title('Gryzinski $\mu_{ion}$ for PMMA components')
plt.xlabel('E, eV')
plt.ylabel('$\mu(E)$, cm$^{-1}$')
plt.legend()
plt.grid()
plt.show()
