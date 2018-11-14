#%% Import
import numpy as np
import importlib
import matplotlib.pyplot as plt
import os
sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
             
os.chdir(sim_path + 'make_e_data')

import sys
sys.path.append(sim_path + 'MODULES')

import my_functions as mf
mf = importlib.reload(mf)
#import matplotlib.pyplot as plt

#%%
def plot_DATA(DATA, d_PMMA=0):
    fig, ax = plt.subplots()
    
    for tn in range(int(np.max(DATA[:, 0]))):
        if len(np.where(DATA[:, 0] == tn)[0]) == 0:
            continue
        beg = np.where(DATA[:, 0] == tn)[0][0]
        end = np.where(DATA[:, 0] == tn)[0][-1] + 1
        
        if np.isnan(DATA[beg, 1]):
            ax.plot(DATA[beg:end, 5], DATA[beg:end, 7], 'k-', linewidth=2)
        else:
            ax.plot(DATA[beg:end, 5], DATA[beg:end, 7], color='0.5',
                    linestyle='dashed', linewidth=2)
        
        ## print events
        inds_el = beg + np.where(DATA[beg:end, 3] == 0)[0]
        inds_ion = beg + np.where(DATA[beg:end, 3] >= 2)[0]
        inds_exc = beg + np.where(DATA[beg:end, 3] == 1)[0]
        ax.plot(DATA[inds_el, 5], DATA[inds_el, 7], 'r*')
        ax.plot(DATA[inds_ion, 5], DATA[inds_ion, 7], 'mo')
        ax.plot(DATA[inds_exc, 5], DATA[inds_exc, 7], 'yv')
    
    points = np.arange(-3e+3, 3e+3, 10)
    ax.plot(points, np.zeros(np.shape(points)), 'k')
    ax.plot(points, np.ones(np.shape(points))*d_PMMA, 'k')
    
    ax.plot(0, 10, 'k-', linewidth=2, label='primary electron')
    ax.plot(0, 10, color='0.5', linestyle='dashed', linewidth=2,
            label='secondary electron')
    ax.plot(0, 10, 'r*', label='elastic scattering')
    ax.plot(0, 10, 'mo', label='ionization')
    ax.plot(0, 10, 'yv', label='excitation')
    
#    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
#    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.gca().set_aspect('equal', adjustable='box')
#    plt.title('Direct Monte-Carlo simulation')
    plt.xlabel('x, nm', fontsize=14)
    plt.ylabel('z, nm', fontsize=14)
    plt.axis('on')
    plt.grid('on')
    plt.xlim(575, 605)
    plt.ylim(970, 1005)
    plt.legend(fontsize=14)
    plt.gca().invert_yaxis()
    
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    
    plt.show()

plot_DATA(DATA, d_PMMA)
