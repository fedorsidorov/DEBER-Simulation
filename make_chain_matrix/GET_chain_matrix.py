#%% Import
import numpy as np
import matplotlib.pyplot as plt
import os
import importlib
import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)
os.chdir(mv.sim_path_MAC + 'make_chain_matrix')

#%%
source_dir = mv.sim_path_MAC + 'CHAINS/950K_122nm/comb_400x100x122/'

source_files = os.listdir(source_dir)

x_min = []
x_max = []

y_min = []
y_max = []

z_min = []
z_max = []

for idx, file in enumerate(source_files):
    
    if file == '.DS_store': continue
    
    mf.upd_progress_bar(idx, len(source_files))
    
    now_chain = np.load(source_dir + file)
    
    x_min.append(now_chain[0, :].min())
    x_max.append(now_chain[0, :].max())
    
    y_min.append(now_chain[1, :].min())
    y_max.append(now_chain[1, :].max())
    
    z_min.append(now_chain[2, :].min())
    z_max.append(now_chain[2, :].max())