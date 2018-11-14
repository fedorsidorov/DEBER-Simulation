#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
os.chdir(sim_path + 'mapping/mass_distr')
import importlib
import my_functions as mf
mf = importlib.reload(mf)
from scipy.interpolate import interp1d

#%%
filenames = ['initial_L', 'final_L_2.5', 'final_L_2.7']
labels = ['initial mass',
          '"2.5" C atoms of 5',
          '"2.7" C atoms of 5']
#colors = ['k', 'deepskyblue', 'darkorange']

final_list = []

for fname in filenames:
    
    final_L = np.load(fname + '.npy')
    for i in range(len(final_L)):
        if final_L[i] == 0:
            final_L[i] = 1
    final_list.append(final_L)

plt.figure()

for i in range(len(final_list)):
    
    log_mw = np.log10(final_list[i] * 100)
    n_bins=25
    mult=1
    if i == 0:
        n_bins = 15
        mult = 0.2
    hist, bins = np.histogram(log_mw, bins=n_bins, normed=True)
    points = (bins[:-1] + bins[1:]) / 2
    plt.plot(points, hist*mult, '-', linewidth=3, label=labels[i])

xrange_B = np.linspace(3.9, 6, num=20)
mat_B = np.loadtxt('curves/sharma_peak_B.dat')
x_log_B = np.log10(mat_B[1:, 0])
y_B = mat_B[1:, 1]
yrange_B = interp1d(x_log_B, y_B)(xrange_B)
plt.bar(xrange_B, yrange_B*3.5 * 0.2,\
        align='center',\
        alpha=0.5,\
        width=0.08,\
        label='peak B')

xrange_A = np.linspace(3.31, 5.27, num=20)
mat_A = np.loadtxt('curves/sharma_peak_A.dat')
x_log_A = np.log10(mat_A[1:, 0])
y_A = mat_A[1:, 1]
yrange_A = interp1d(x_log_A, y_A)(xrange_A)
plt.bar(xrange_A, yrange_A*0.85,\
        align='center',\
        alpha=0.5,\
        width=0.08,\
        label='peak A')

plt.xlim([3, 6])
plt.ylim([0, 1])
plt.legend(loc='upper left')
plt.xlabel('log(M$_W$)')
plt.ylabel('Arbitrary units')
plt.title('Simulation of mass distribution')
plt.grid()
plt.show()
