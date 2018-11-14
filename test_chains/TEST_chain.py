#%% Import
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append('../MODULES')

import importlib

import my_arrays_x_smart as ma
import my_variables as mv
import my_functions as mf

ma = importlib.reload(ma)
mv = importlib.reload(mv)
mf = importlib.reload(mf)

#%%
my_bins = np.array((1, 2, 3, 4, 5, 6))

#arr += np.histogramdd(10, bins=my_bins)[0]

chain_arr = np.ones((10, 3))*1

cube_5x5 = np.zeros((5, 5, 5))

pos = chain_arr[0, :]

cube_5x5 += np.histogramdd(chain_arr, bins=(my_bins, my_bins, my_bins))[0]

#%%
def filter_chain_with_5x5_cube(chain_arr, xyz):
    
    step = 0.25e-7
    my_bins = []
    
    ## cell edges of 5x5 cube
    for coord in xyz:
        my_bins.append(np.arange(-2.5, 2.5 + 0.1)*step + coord)

    ## 5x5 cube histogram
    cube_5x5 = np.zeros((5, 5, 5))
    
    ## check all the links from 0 to l for current chain
    for i in range(len(chain_arr)):
        pos = chain_arr[i, :] ## current link position
        
        if np.all([my_bins[i].min() < pos[i] < my_bins[i].max() for i in range(3)]):
            cube_5x5 += np.histogramdd(pos.reshape((1, 3)), bins=my_bins)[0]
                            
    return cube_5x5

def get_cube_5x5_list(cube_5x5, xyz):
    
    step = 0.25e-7
    my_centers = []
    
    ## cell centers of 5x5 cube
    for coord in xyz:
        my_centers.append(np.arange(-2, 2 + 0.1)*step + coord)
        
    ## list with possible positions for the next link
    cube_5x5_list = []
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                if i == j == k == 2:
                    continue
                ## if there is only one neighbour (previous link), OK
                if np.sum(cube_5x5[i-1:i+2, j-1:j+2, k-1:k+2]) == 1:
                    cube_5x5_list.append((my_centers[0][i], my_centers[1][j], my_centers[2][k]))
    
    return cube_5x5_list

def make_chain(start_pos, L, chains_list):
    
    print('Make chain with L =', L)
    chain_arr = np.zeros((L, 3))*np.nan
    
    ## chain starting point
    xyz = start_pos
    chain_arr[0, :] = xyz
    
    ## collision link number
    jam_pos = 0
    ## collision counter
    jam_cnt = 0
    
    l = 1
    
    while l < L:
        
        mf.upd_progress_bar(l, L)
        
        ## get 5x5 cube with existent monomers
        cube_5x5 = np.zeros((5, 5, 5))
        
        for chain in [*chains_list, chain_arr]:
            cube_5x5 += filter_chain_with_5x5_cube(chain, xyz)
        
        ## get list of 5x5 cube with appropriate positions for the next monomer
        cube_5x5_list = get_cube_5x5_list(cube_5x5, xyz)

        ## if no free space
        if len(cube_5x5_list) == 0:
                        
            if np.abs(jam_pos - l) < 10: ## if new jam is near current link
                jam_cnt += 1 ## increase collision counter
            
            else: ## if new jam is on new link
                jam_pos = l ## set new collision position
                jam_cnt = 0 ## set collision counter to 0
            
            print(l, ': No free space,', jam_cnt)
            
            ## if possible, make rollback proportional to jam_cnt
            rollback_step = jam_cnt // 5 + 1
            
            if l - (rollback_step + 1) >= 0:
                xyz = chain_arr[l - (rollback_step + 1), :]
                l -= rollback_step
                continue
            
            else:
                print('Jam in very start!')
                xyz = start_pos
                l = 1
                continue

        ## define pos in cube_5x5_list for next link
        new_link = cube_5x5_list[mf.randint(len(cube_5x5_list))]
        chain_arr[l, :] = new_link
        xyz = new_link
        
        l += 1
    
    chains_list.append(chain_arr)
    
    return chain_arr, chains_list

#%%
resist_arr = np.ones((400, 400, 640))*(-1)
L_arr = np.load('L_arr_7_jul.npy')

chains_list = []

## сhain step
step = 0.25e-7

## chain_num
chain_num = 0

for i in range(400):
    for j in range(400):
        for k in range(640):
            
            L = int(L_arr[chain_num % len(L_arr)])
            start_pos = np.array((i, j, k))*step
            chain_arr, chains_list = make_chain(start_pos, L, chains_list)
            
            chain_num += 1
            
#%%
mf.print_chains_list(chains_list)
