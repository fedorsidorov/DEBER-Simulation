#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
#import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'test_chain_size')

#%%
source_dir = mv.sim_path_MAC + 'CHAINS/160_mon/'

#N_chains = 10000
#
#chain_bank = []
#
#for i in range(N_chains):
#    
#    mf.upd_progress_bar(i, N_chains)
#    
#    now_chain = np.load(source_dir + 'chain_' + str(i) + '.npy')
#    chain_bank.append(now_chain)

#%%
n_mon = 160
l_x = l_y = l_z = 100e-7

vol = l_x * l_y * l_z

n_mon_required = vol * mv.rho_PMMA / mv.m_PMMA_mon
n_chains_required = n_mon_required / n_mon

chain_size = (vol / n_chains_required) ** (1/3) / 1e-7 ## 0.28 nm


