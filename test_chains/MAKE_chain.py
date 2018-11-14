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
mat = np.loadtxt('curves/sharma_peak_B.dat')
x_log = np.log10(mat[1:, 0])
y = mat[1:, 1]
X_LOG = np.arange(x_log[0], x_log[-1], 1e-2)
Y_LOG = mf.log_interp1d(x_log, y)(X_LOG)
#plt.plot(X_LOG, Y_LOG)
plt.semilogx(X_LOG, Y_LOG, 'ro')
plt.show()

#%%
S_arr = np.zeros(np.size(X_LOG))
s_tot = np.trapz(Y_LOG, x=X_LOG)

for i in range(len(S_arr) - 1):
    S_arr[i] = np.trapz(Y_LOG[0:i+1], x=X_LOG[0:i+1])/s_tot

S_arr[-1] = 1

plt.plot(X_LOG, S_arr)
plt.show()

#%%
def get_log_mw():
    r = mf.random()
    for i in range(len(S_arr) - 1):
        if r < S_arr[i + 1]:
            return X_LOG[i]

def get_z0():
    return 0
#    r = mf.random()
#    return (10 + r*1180)*1e-7 - 600e-7
    #return 600e-7

#%%
N_chains = 10000

log_mw_arr = np.zeros(N_chains)
L_arr = np.zeros(N_chains)

for i in range(N_chains):
    log_mw = get_log_mw()
    log_mw_arr[i] = log_mw
    L_arr[i] = int(10**log_mw / 100)

z0_arr = np.zeros(N_chains)
for i in range(N_chains):
    z0 = get_z0()
    z0_arr[i] = (z0 // 0.25e-7) * 0.25e-7

#%%
L_arr = np.load("./L_arr_7_jul.npy")
mf.print_histogram(np.log10(L_arr*100), 30)
plt.semilogx(X_LOG, Y_LOG, 'ro')

#%%
## сhain step
step = 0.25e-7

## chain_num
#chain_num = 0
chain_num = 7006

while chain_num < N_chains:
    
    ## current chain length
    L = int(L_arr[chain_num])
    chain_arr = np.zeros((L, 3))*np.nan
    
#    z0 = z0_cube_5x5[chain_num]
    z0 = 0
    (x, y, z) = (0, 0, z0)
    
    print('L =', L, 'z0 =', int(z0*1e+7))
    
    ## chain starting point
    chain_arr[0, :] = (x, y, z)
    l = 1
    
    ## collision counter
    jam_cnt = 0
    ## collision link number
    jam_pos = 0
    
    while l < L:      

        ## cell edges of 5x5 cube
        x_e = np.arange(x - step*2.5, x + step*2.5 + 0.01*step, step)
        y_e = np.arange(y - step*2.5, y + step*2.5 + 0.01*step, step)
        z_e = np.arange(z - step*2.5, z + step*2.5 + 0.01*step, step)
        
        ## cell centers of 5x5 cube
        x_c = np.arange(x - step*2, x + step*2 + 0.01*step, step)
        y_c = np.arange(y - step*2, y + step*2 + 0.01*step, step)
        z_c = np.arange(z - step*2, z + step*2 + 0.01*step, step)
        
        ## 5x5 cube histogram
        cube_5x5 = np.zeros((5, 5, 5))
        
        ## check all the links from 0 to l
        for i in range(l):
            p = chain_arr[i, :] # current link position
            if x_e[0] < p[0] < x_e[-1]:
                if y_e[0] < p[1] < y_e[-1]:
                    if z_e[0] < p[2] < z_e[-1]:
                        ## if any link is in 5x5 cube, add it to cube_5x5
                        cube_5x5 += np.histogramdd(p.reshape((1, 3)), \
                                    bins=(x_e, y_e, z_e))[0]
        
        ## list with possible positions for the next link
        cube_5x5_list = []
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    if i == j == k == 2:
                        continue
                    ## if z is out of range, continue
#                    if z_c[k] < 0:
#                        continue
                    ## if there is only one neighbour (previous link), OK
                    if np.sum(cube_5x5[i-1:i+2, j-1:j+2, k-1:k+2]) == 1:
                        cube_5x5_list.append((x_c[i], y_c[j], z_c[k]))
        
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
                (x, y, z) = chain_arr[l - (rollback_step + 1), :]
                l -= rollback_step
                continue
            
#            if l - (2 + jam_cnt // 5) >= 0:
#                (x, y, z) = chain_arr[l - (2 + jam_cnt // 5), :]
#                l -= (1 + jam_cnt // 5)
#                continue
            
            else:
                print('Jam in very start!')
                break
        
        ## define pos in cube_5x5_list for next link
#        pos = mf.randint(0, len(cube_5x5_list) - 1)
#        new_link = cube_5x5_list[pos]
        new_link = mf.choice(cube_5x5_list)
        chain_arr[l, :] = new_link
        (x, y, z) = new_link
        
        mf.upd_progress_bar(L, l + 1)
        l += 1
     
    if l < L:
        continue
    
    else: ## get some info about chain
#        x_min_str = str(int(np.min(chain_arr[:, 0])*1e+7))
#        x_max_str = str(int(np.max(chain_arr[:, 0])*1e+7))
#        y_min_str = str(int(np.min(chain_arr[:, 1])*1e+7))
#        y_max_str = str(int(np.max(chain_arr[:, 1])*1e+7))
#        z_min_str = str(int(np.min(chain_arr[:, 2])*1e+7))
#        z_max_str = str(int(np.max(chain_arr[:, 2])*1e+7))
        
        fname_list = ['HOME_chains_7_jul/Ch',\
                      str(chain_num)]
#                      'z0=' + str(int(z0*1e+7)), 'chain_num=' + \
#            str(chain_num), 'L=' + str(L), \
#            'x', x_min_str, x_max_str, 'y', y_min_str, y_max_str, \
#            'z', z_min_str, z_max_str]
        fname = '_'.join(fname_list)
        
        ## convert cm to nm
        chain_arr_nm = chain_arr*1e+7
        result = np.array(np.hstack((chain_arr_nm, np.array(range(L)).reshape((L, 1)))))
        
        np.save(fname, result)
        print('Chain ' + str(chain_num) + ' is saved, file:', fname)  
        chain_num += 1

#%%
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
chain = np.load('HOME_chains_7_jul/Ch_3534.npy')
#ax.plot(chain[:, 0], chain[:, 1], chain[:, 2])
ch_x = chain[:, 0]

#%%
#beg = 0
#end = 4572
beg = 0
end = 15493

ax.plot(chain_arr_nm[beg:end, 0], chain_arr_nm[beg:end, 1], chain_arr_nm[beg:end, 2], label='PMMA chain, L = 15000')

ax.plot(chain_arr_nm[99:100, 0], chain_arr_nm[99:100, 1], chain_arr_nm[99:100, 2], 'ro', label='scissions')
ax.plot(chain_arr_nm[100:1100, 0], chain_arr_nm[100:1100, 1], chain_arr_nm[100:1100, 2], color='orange', label='detached monomers')

ax.plot(chain_arr_nm[599:600, 0], chain_arr_nm[599:600, 1], chain_arr_nm[599:600, 2], 'ro')
ax.plot(chain_arr_nm[600:1600, 0], chain_arr_nm[600:1600, 1], chain_arr_nm[600:1600, 2], color='orange')

ax.plot(chain_arr_nm[3099:3100, 0], chain_arr_nm[3099:3100, 1], chain_arr_nm[3099:3100, 2], 'ro')
ax.plot(chain_arr_nm[3100:4100, 0], chain_arr_nm[3100:4100, 1], chain_arr_nm[3100:4100, 2], color='orange')

ax.plot(chain_arr_nm[5499:5500, 0], chain_arr_nm[5499:5500, 1], chain_arr_nm[5499:5500, 2], 'ro')
ax.plot(chain_arr_nm[5500:6500, 0], chain_arr_nm[5500:6500, 1], chain_arr_nm[5500:6500, 2], color='orange')

ax.plot(chain_arr_nm[8099:8100, 0], chain_arr_nm[8099:8100, 1], chain_arr_nm[8099:8100, 2], 'ro')
ax.plot(chain_arr_nm[8100:9100, 0], chain_arr_nm[8100:9100, 1], chain_arr_nm[8100:9100, 2], color='orange')

ax.plot(chain_arr_nm[9099:9100, 0], chain_arr_nm[9099:9100, 1], chain_arr_nm[9099:9100, 2], 'ro')
ax.plot(chain_arr_nm[9100:10100, 0], chain_arr_nm[9100:10100, 1], chain_arr_nm[9100:10100, 2], color='orange')

ax.plot(chain_arr_nm[13199:13200, 0], chain_arr_nm[13199:13200, 1], chain_arr_nm[13199:13200, 2], 'ro')
ax.plot(chain_arr_nm[13200:14200, 0], chain_arr_nm[13200:14200, 1], chain_arr_nm[13200:14200, 2], color='orange')

ax.set_xlabel('x, nm')
ax.set_ylabel('y, nm')
ax.set_zlabel('z, nm')
plt.legend()
plt.show()

#%%
res_cube_5x5 = np.load('Ch_z0=0_chain_num=27_L=2232_x_-9_17_y_-19_12_z_0_16.npy')



dtype = [('x', float), ('y', float), ('z', float), ('n', int)]
res_cube_5x5_v = np.cube_5x5ay(res_cube_5x5, dtype=dtype)
res_cube_5x5_v_sort = np.sort(res_cube_5x5_v, axis=0, order=['x', 'y', 'z'])

#%%
L_cube_5x5_slice = L_cube_5x5[769:3483]
M_slice = np.log10(L_cube_5x5_slice * 100)

