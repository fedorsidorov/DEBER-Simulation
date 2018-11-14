import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit
import sys
sys.path.append('../MODULES')
import my_functions as mf

#%% Load
eps_real_raw = np.loadtxt('curves/eps_real.txt')
eps_imag_raw = np.loadtxt('curves/eps_imag.txt')

plt.plot(eps_real_raw[:, 0], eps_real_raw[:, 1], 'r',\
         eps_imag_raw[:, 0], eps_imag_raw[:, 1], 'b--')

hw = np.arange(0, 40, 0.1)

#%% Interpolate
lin_interp_real = interp1d(eps_real_raw[:, 0], eps_real_raw[:, 1], kind='linear')
lin_interp_imag = interp1d(eps_imag_raw[:, 0], eps_imag_raw[:, 1], kind='linear')

eps_real = lin_interp_real(hw)
eps_imag = lin_interp_imag(hw)

plt.plot(hw, eps_real, 'r', hw, eps_imag, 'b--')

#%%
#def get_Im(re, im):
#    return im / (re**2 + im**2)
#
#plt.plot(hw, get_Im(eps_real, eps_imag), 'r')

#%% Get ELF
Im_raw = np.loadtxt('curves/Im.txt')
plt.plot(Im_raw[:, 0], Im_raw[:, 1], 'b')
#
#def Im_tail_func(x, a, b, c, d):
#        return  a + b / x + c * np.exp(-d*x)
#
#def fit_Im_tail():
#
#    Im_raw = np.loadtxt('curves/Im.txt')
#    start_ind = 70
#    popt, pcov = curve_fit(Im_tail_func, Im_raw[start_ind:, 0], Im_raw[start_ind:, 1])    
#
#    return popt
#
##hw_ext = np.arange(0, 100, 0.1)
#hw_ext = np.arange(0, 15000, 0.5)
#
#plt.plot(hw_ext, Im_tail_func(hw_ext, *fit_Im_tail()), 'r-')
#plt.show()

#%%
#hw_ext = np.arange(0, 100, 0.1)
hw_ext = np.arange(0, 15000, 0.5)

lin_interp_Im_raw = interp1d(Im_raw[:, 0], Im_raw[:, 1], kind='linear')
interm_ind = 120

Im_ext = np.zeros(len(hw_ext))

Im_ext[:interm_ind] = lin_interp_Im_raw(hw_ext[:interm_ind])
#Im_ext_end = Im_tail_func(hw_ext[interm_ind:], *fit_Im_tail())

#Im_ext = np.hstack((Im_ext_beg, Im_ext_end))

#Im_ext_filt = medfilt(Im_ext, kernel_size=35)


#plt.plot(hw_ext[:1000], Im_ext[:1000], 'b')
#plt.plot(hw_ext[:1000], Im_ext_filt[:1000], 'g')
#plt.show()

#%% Plan B
def get_Im_theor(hw):
    
#    h = 4.13e-15
    
    A = np.array((129.22, 60.45, 159.60))
    hw0 = np.array((19.7, 23.8, 30))
    g = np.array((10.5, 7.2, 13))
    Eg = 5
#    g = hg / h
    
    beg_ind = np.where(hw > Eg)[0][0]
    
    Im_theor = np.zeros(len(hw))
    
    for i in range(3):
        
        Im_theor[beg_ind:] += A[i]*g[i]*hw[beg_ind:] /\
            ((hw0[i] - hw[beg_ind:]**2)**2 + (g[i]*hw[beg_ind:])**2)
        
    return Im_theor

Im_theor = get_Im_theor(hw_ext)

#plt.plot(hw_ext[:1000], Im_theor[:1000])

Im_ext[interm_ind:] = Im_theor[interm_ind:]#*0.96/6.42

plt.plot(Im_raw[:, 0], Im_raw[:, 1], 'ro')
plt.plot(hw_ext[:1500], Im_ext[:1500], 'g')
plt.show

#%% Get dEds
def v(a):
    s = np.sqrt(1 - 2*a)
    A = 2*s / ((1 + a)*(1 + a + s))
    B = ((1 - a**2)*(1 + a)) / ((1 - a - s)*(1 + a + s)**2)
    return A + np.log(B)

def w(a):
    s = np.sqrt(1 - 2*a)
    A = (3 * a**2 + 3*a + 1) / ((1 + a)**2)
    B = (1 + a - s) / (1 + a)
    C = (1 - a) / (1 - a - s)
    D = (2 * a**2 + a) / ((1 + a)**2)
    E = (1 + a) / (1 + a + s)
    F = (2*a*s) / ((1 + a)**2 * (1 + a + s))
    return A*np.log(B) + np.log(C) + D*np.log(E) + F

a0 = 5.3e-9

E_arr = np.arange(0, 30000, 0.5)
SP = np.zeros(len(E_arr))
IMFP_inv = np.zeros(len(E_arr))

for i in range(len(E_arr)):
    
    mf.upd_progress_bar(i, len(E_arr))
    
    E = E_arr[i]
    
    if E == 0:
        SP[i] = 0
        IMFP_inv[i] = 1e+10
        continue
    
    inds = np.where(E_arr <= E/2)[0][1:]
    hw_arr = E_arr[inds]
    
    Im_arr = Im_ext[inds]
    alpha_arr = hw_arr / E
    v_arr = v(alpha_arr)
    w_arr = w(alpha_arr)
    
    SP_integral_arr = hw_arr * Im_arr * v_arr
    IMFP_inv_integral_arr = Im_arr * w_arr
    
    SP[i] = 1/(2 * np.pi * a0 * E) * np.trapz(SP_integral_arr, x=hw_arr, dx=0.5)
    IMFP_inv[i] = 1/(2 * np.pi * a0 * E) * np.trapz(IMFP_inv_integral_arr, x=hw_arr, dx=0.5)
    
#%%
SP_paper = np.loadtxt('curves/SP.txt')

plt.figure(1)

plt.loglog(SP_paper[:, 0], SP_paper[:, 1], label='paper')
plt.loglog(E_arr, SP*1e-8, label='me')
plt.xlim(2e+2, 3e+4)
plt.ylim(3e-2, 2)
plt.xlabel('eV')
plt.ylabel('dE/ds, ev/A')
plt.legend()
plt.show()

#%%
IMFP_paper = np.loadtxt('curves/IMFP.txt')

plt.figure(2)

plt.loglog(IMFP_paper[:, 0], IMFP_paper[:, 1], label='paper')
plt.loglog(E_arr, 1/IMFP_inv*1e+8, label='me')
#plt.xlim(2e+2, 3e+4)
#plt.ylim(3e-2, 2)
plt.xlabel('eV')
plt.ylabel('IMFP, A')
plt.legend()
plt.show()

