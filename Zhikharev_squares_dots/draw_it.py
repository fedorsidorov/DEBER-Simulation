#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'Zhikharev_squares_dots')

#%%
dots = np.loadtxt('dots_cross_1.txt')

plt.figure(figsize=[6, 2.2])

plt.plot(dots[:, 0]*1e+7, dots[:, 1]*1e+7, label='scan 1')
plt.plot(dots[:, 2]*1e+7, dots[:, 3]*1e+7, label='scan 2')

plt.ylim(0, 2.5)
#plt.ylim(0.75, 2.25)
#plt.legend(loc='4', fontsize=14)

ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)

plt.xlabel('x, nm', fontsize=14)
plt.ylabel('z, nm', fontsize=14)

ax.tick_params(axis='x',direction='in', pad=-15)

plt.grid()

#%%
squares = np.loadtxt('squares_cross_2.txt')

plt.figure(figsize=[6, 2.2])

plt.plot(squares[:, 0]*1e+7, squares[:, 1]*1e+7, label='scan 1')
plt.plot(squares[:, 2]*1e+7, squares[:, 3]*1e+7, label='scan 2')
plt.plot(squares[:, 4]*1e+7, squares[:, 5]*1e+7, label='scan 3')
plt.plot(squares[:, 6]*1e+7, squares[:, 7]*1e+7, label='scan 4')

plt.xlim(0, 175)
plt.ylim(0, 3.5)
#plt.legend(loc='4', fontsize=14)

ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)

plt.xlabel('x, nm', fontsize=14)
plt.ylabel('z, nm', fontsize=14)

ax.tick_params(axis='x',direction='in', pad=-15)

plt.grid()
