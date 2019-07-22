#%% Import
import numpy as np
import os
import importlib
import matplotlib.pyplot as plt

import my_functions as mf
import my_variables as mv
import SAW_chain as saw

mf = importlib.reload(mf)
mv = importlib.reload(mv)
saw = importlib.reload(saw)

os.chdir(mv.sim_path_MAC + 'test_chain_size')

#%%
n_chains = 10000
chain_num = 0

chain_list = []

while chain_num < n_chains:
    
    chain_len = 160
    
    now_chain = saw.make_PMMA_chain(chain_len)
    
    dirname = '160_mon/'
    filename = 'chain_' + str(chain_num) + '.npy'
    
    np.save(dirname + filename, now_chain)
    print(filename + ' is saved')
    
    chain_list.append(now_chain)
    chain_num += 1
