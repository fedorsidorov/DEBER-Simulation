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

os.chdir(mv.sim_path_MAC + 'diel_responce')

#%%
E_arr = np.logspace(0, 4.4, 10000)
#E_arr = np.linspace(1, 20e+3, 10000)

Im_arr = np.zeros(len(E_arr))

## En, Gn, An from dapor2015.pdf
params = [
        [19.46, 8.770, 100.0],
        [25.84, 14.75, 286.5],
        [300.0, 140.0, 80.0],
        [550.0, 300.0, 55.0],
        ]

for arr in params:
    
    E, G, A, = arr
    Im_arr += A*G*E_arr / ((E**2 - E_arr**2)**2 + (G*E_arr)**2)

plt.loglog(E_arr, Im_arr, 'r.', label='oscillators')

plt.legend()
plt.xlabel('E, eV')
plt.ylabel('Im[1/eps]')
plt.grid()
plt.show()

#%%
def L(x):
    f = (1-x)*np.log(4/x) - 7/4*x + x**(3/2) - 33/32*x**2
    return f

def S(x):
    f = np.log(1.166/x) - 3/4*x - x/4*np.log(4/x) + 1/2*x**(3/2) - x**2/16*np.log(4/x) - 31/48*x**2
    return f

#%%
IMFP_inv_arr = np.zeros(len(E_arr))
dEds_arr = np.zeros(len(E_arr))

diff_IMFP_inv_arr = np.zeros((len(E_arr), len(E_arr)))

for i in range(len(E_arr)):
    
    now_E = E_arr[i]
    
    inds = np.where(E_arr <= now_E/2)[0]
    
    y_IMFP = Im_arr[inds] * L(E_arr[inds]/E_arr[i])
    y_dEds = Im_arr[inds] * S(E_arr[inds]/E_arr[i]) * E_arr[inds]*mc.eV
    
    IMFP_inv_arr[i] = mc.k_el * mc.m * mc.e**2 / (2 * np.pi * mc.hbar**2 * E_arr[i]*mc.eV) *\
        np.trapz(y_IMFP, x=E_arr[inds]*mc.eV)
        
    dEds_arr[i] = mc.k_el * mc.m * mc.e**2 / (np.pi * mc.hbar**2 * E_arr[i]*mc.eV) *\
        np.trapz(y_dEds, x=E_arr[inds]*mc.eV)
    
    diff_IMFP_inv_arr[i, inds] = mc.k_el * mc.m * mc.e**2 / (2 * np.pi * mc.hbar**2 *\
        E_arr[i]*mc.eV) * Im_arr[inds] * L(E_arr[inds]/E_arr[i])
    

#%%
IMFP_solid = np.loadtxt('curves/IMFP_solid.txt')
IMFP_dashed = np.loadtxt('curves/IMFP_dashed.txt')

plt.loglog(IMFP_solid[:, 0], IMFP_solid[:, 1], label='Dapor_solid')
plt.loglog(IMFP_dashed[:, 0], IMFP_dashed[:, 1], label='Dapor_dashed')

plt.loglog(E_arr, 1/IMFP_inv_arr * 1e+10, label='My')

plt.xlim(10, 10000)
plt.ylim(1, 1000)

plt.xlabel('E, eV')
plt.ylabel('IMFP, $\AA$')
plt.legend()
plt.grid()
plt.show()

#%%
dEds_solid = np.loadtxt('curves/dEds_solid.txt')
dEds_dashed = np.loadtxt('curves/dEds_dashed.txt')
dEds_dotted = np.loadtxt('curves/dEds_dotted.txt')

plt.semilogx(dEds_solid[:, 0], dEds_solid[:, 1], label='Dapor_solid')
plt.semilogx(dEds_dashed[:, 0], dEds_dashed[:, 1], label='Dapor_dashed')
plt.semilogx(dEds_dotted[:, 0], dEds_dotted[:, 1], label='Dapor_dotted')

plt.semilogx(E_arr, dEds_arr / mc.eV / 1e+10, label='My')
plt.xlim(10, 10000)
plt.ylim(0, 4)

plt.xlabel('E, eV')
plt.ylabel('dEds, eV/$\AA$')
plt.legend()
plt.grid()
plt.show()

#%%
for i in range(0, len(E_arr), 1000):
    
    plt.loglog(E_arr, diff_IMFP_inv_arr[i, :])

#%%
IMFP_inv_arr_test = np.zeros(len(E_arr))

for i in range(len(E_arr)):
    
    now_E = E_arr[i]
    
    inds = np.where(E_arr <= now_E/2)[0]
    
    y_IMFP = Im_arr[inds] * L(E_arr[inds]/E_arr[i])
    y_dEds = Im_arr[inds] * S(E_arr[inds]/E_arr[i]) * E_arr[inds]*mc.eV
    
    IMFP_inv_arr_test[i] = np.trapz(diff_IMFP_inv_arr[i, :], x=E_arr*mc.eV)

#%%
plt.loglog(E_arr, 1/IMFP_inv_arr_test * 1e+10, label='My')
