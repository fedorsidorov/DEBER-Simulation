#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os

sim_path = '/Users/fedor/.yandex.disk/434410540/Yandex.Disk.localized/' +\
            'Study/Simulation/'
#sim_path = '/home/fedor/Yandex.Disk/Study/Simulation/'
             
os.chdir(sim_path + 'mapping')

#import sys
#sys.path.append(sim_path + 'MODULES')

import importlib

import my_functions as mf
mf = importlib.reload(mf)

#%%
def get_x0y0(lx, ly, space):
    return mf.uniform(-space, lx + space), mf.uniform(-space, ly + space)


source_folder = '../make_e_data/DATA_5keV_600nm/'
e_data_filenames = os.listdir(source_folder)
dest_folder = 'Courtney/DATA_Courtney/'

n_files = 0
n_files_max = len(e_data_filenames)
#n_files_max = 500

lx = 100
ly = 100
space = 2

#%%
for fname in e_data_filenames:
    
    mf.upd_progress_bar(n_files, n_files_max)
    
    if fname == '.DS_Store':
        continue
    
    DATA_Pn = np.load(source_folder + fname)
    
    for i in range(5):
    
        now_DATA_Pn = DATA_Pn.copy()
        
        ## make rotation
        now_DATA_Pn[:, 5:7] = mf.add_xy_rotation(now_DATA_Pn[:, 5:7], 2*np.pi*mf.random())
        
        ## primary tracks number
        n_tr_prim = int(DATA_Pn[np.where(np.isnan(DATA_Pn[:, 1]))][-1, 0] + 1)
        
        for track_num in range(n_tr_prim):
        
            ## in case of only elastic events in PMMA
            if len(np.where(DATA_Pn[:, 0] == track_num)[0]) == 0:
                print('No events')
                print(fname, track_num)
                continue
            
            ## in normal case
            else:
                x0, y0 = get_x0y0(lx, ly, space)            
                now_DATA_Pn = mf.add_xy_shift(now_DATA_Pn, track_num, x0, y0)
        
        np.save(dest_folder + fname.replace('.npy', '_Courtney_' + str(i) + '.npy'), now_DATA_Pn)

    n_files += 1

#%% check
x_list = []
y_list = []

for fname in os.listdir(dest_folder):
    
    if fname == '.DS_Store':
        continue
    
    DATA_Pn = np.load(dest_folder + fname)
    
    for line in DATA_Pn:
        if line[7] == 0 and line[4] == 25000:
            
            x_list.append(line[5])
            y_list.append(line[6])
