#%% Import
import numpy as np
import os
import importlib
import my_functions as mf
import my_variables as mv
import matplotlib.pyplot as plt

mf = importlib.reload(mf)
mv = importlib.reload(mv)

os.chdir(mv.sim_path_MAC + 'test_chain_size')

#%%
def check_chain(chain_coords, now_mon_coords, d_2):
    for mon_coords in chain_coords[:-1, :]:    
        if np.sum((mon_coords - now_mon_coords)**2) < d_2:
            return False
    return True

#%%
d_mon = 0.28
d_mon_2 = d_mon**2

#lz = 500

theta = np.deg2rad(109)

#n_chains = len(L_arr)
#n_chains = 1000
n_chains = 1

chain_num = 0

chains_list = []

while chain_num < n_chains:
    
    L = 160
    print('New chain, L =', L)
    
    chain_coords = np.zeros((L, 3))
    
    chain_coords[0, :] = 0, 0, 0
    
    On = mf.get_O_matrix(2 * np.pi * mf.random(), theta, np.eye(3))
    
    ## collision counter
    jam_cnt = 0
    ## collision link number
    jam_pos = 0
    
    i = 1
    
    while i < L:
        
        mf.upd_progress_bar(i, L)
            
        while True:
            
            dxdydz = On.transpose() * np.mat([[0], [0], [1]]) * d_mon
            chain_coords[i, :] = chain_coords[i-1, :] + dxdydz.A1
            On = mf.get_O_matrix(2 * np.pi * mf.random(), theta, On)
            
            st = check_chain(chain_coords[:i, :], chain_coords[i, :], d_mon_2)
                
            if st:
                break
            
            else: ## if no free space
                
                if np.abs(jam_pos - i) < 10: ## if new jam is near current link
                    jam_cnt += 1 ## increase collision counter
                
                else: ## if new jam is on new link
                    jam_pos = i ## set new collision position
                    jam_cnt = 0 ## set collision counter to 0
                
                print(i, ': No free space,', jam_cnt)
                
                ## if possible, make rollback proportional to jam_cnt
                rollback_step = jam_cnt // 10
                
                if i - (rollback_step + 1) >= 0:
                    i -= rollback_step
                    continue
                
                else:
                    print('Jam in very start!')
                    break
        
        i += 1
    
#    dirname = '160_mon/'
#    filename = 'chain_' + str(chain_num) + '.npy'
#    np.save(dirname + filename, chain_coords)
#    print(filename + ' is saved')
#    chains_list.append(chain_coords)
    chain_num += 1
