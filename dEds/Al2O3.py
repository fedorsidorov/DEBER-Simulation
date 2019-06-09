#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'dEds')

#%%
rho = 3.97e+3 ## kg / m^3
eV = 1.6e-19 ## J

def get_Bethe_dEds(E):
    
    ## Al2O3

    Z = 50
#    A = 102
    Na = 6.02e+23
    I = 145.2 ## eV
#    I = 515 ## eV
    e = 1.6e-19
    k = 9e+9
    eV = 1.6e-19 ## J
    u = 102e-3 ## kg / mole
    Z = 50
    
#    ZA = 0.49    
#    mc2 = 511e+3 ## eV
#    t = E / mc2
#    delta = 2 * np.log(28.816 * np.sqrt(rho * ZA)) + 2 * np.log(t + 1)
#    B0 = np.log(t**2 * (t + 2) / 2) + (1 + t**2 / 8 - (2 * t + 1) * np.log(2)) / (t + 1)**2
#    B = B0 - 2 * np.log(I / mc2) - delta
#    beta2 = 2 * E / mc2
#    dEds = rho * 0.153536 / beta2 * ZA * B
    
    dEds = k**2 * 2 * np.pi * e**4 * (rho * Na / u) * Z * (1 / (E * eV)) * np.log(1.66 * E / I) ## J
    
    return dEds

#%%
Akk_arr = np.loadtxt('Al2O3/Akkerman.txt')
Im_arr = np.loadtxt('Al2O3/Im_eps_corr.txt')

plt.loglog(Akk_arr[:, 0], Akk_arr[:, 1], 'ro')
plt.loglog(Im_arr[:, 0], Im_arr[:, 1])

E_Bethe = np.power(10, np.linspace(1, 4, 100))
dEds_Bethe = get_Bethe_dEds(E_Bethe) * 1e+8 ##


plt.loglog(E_Bethe, get_Bethe_dEds(E_Bethe) / eV / 1e+10, label='Bethe')

plt.xlabel('E, эВ')
plt.ylabel(r'$\frac{dE}{ds}$, эВ / $\AA$')

plt.xlim(10, 4e+4)
plt.ylim(0.1, 10)

#plt.legend()

plt.grid()
plt.show()

#plt.savefig('dEds.png', dpi=300)
