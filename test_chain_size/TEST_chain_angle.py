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

import math

#%%
def check_chain(chain_coords, now_mon_coords, d_2):
    for mon_coords in chain_coords[:-1, :]:    
        if np.sum((mon_coords - now_mon_coords)**2) < d_2:
            return False
    return True

#%%
def get_d_xyz(L, PT, PT_old):
    
    P, T = PT[0], PT[1]
    P_old, T_old = PT_old[0], PT_old[1]
    
    if T_old == 0:
        T_old = 1e-5
    
    if T_old == np.pi:
        T_old = np.pi - 1e-5
    
    cosT_new = np.cos(T_old)*np.cos(T) + np.sin(T_old)*np.sin(T)*np.cos(P) 
    sinT_new = np.sin(np.arccos(cosT_new))
    
    print(sinT_new)
    print(cosT_new)
    
    T_new = np.arccos(cosT_new)
    
    if sinT_new == 0:
        sinT_new = 1e-5
    
    cosDP = (np.cos(T) - cosT_new*np.cos(T_old)) / (np.sin(T_old)*sinT_new)
    sinDP = np.sin(T)*np.sin(P) / sinT_new
    
    a = np.array([[np.cos(P_old), -np.sin(P_old)], [np.sin(P_old), np.cos(P_old)]])
    b = np.array([cosDP, sinDP])
    
    sinP_new, cosP_new = np.linalg.solve(a, b)
    
    print(cosP_new)
    
    if cosP_new > 1:
        cosP_new = 1
    
    if cosP_new < -1:
        cosP_new = -1
    
    P_new = np.arccos(cosP_new)
    
    if sinP_new < 0:
        P_new = -P_new
    
    dx = sinT_new*cosP_new
    dy = sinT_new*sinP_new
    dz = L*cosP_new
    
    return dx, dy, dz, P_new, T_new
    
#%%
d_mon = 0.28
d_mon_2 = d_mon**2

theta_PMMA = np.deg2rad(109)

#n_chains = 1000
n_chains = 1

chain_num = 0

chains_list = []

chain_arr = []

while chain_num < n_chains:
    
    L = 160
    
    print('New chain, L =', L)
    
    chain_coords = np.zeros((L, 3))
    
    chain_coords[0, :] = 0, 0, 0
    
    phi, theta = 0, 0
    
    ## collision counter
    jam_cnt = 0
    ## collision link number
    jam_pos = 0
    
    i = 1
    
    while i < L:
        
        mf.upd_progress_bar(i, L)
        
        while True:
            
            P, T = 0, theta_PMMA
            
            dx, dy, dz, phi, theta = get_d_xyz(d_mon, [P, T], [phi, theta])
            chain_coords[i, :] = chain_coords[i-1, :] + np.array((dx, dy, dz))
            
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
    
    chain_arr = chain_coords
    
    chain_num += 1

#%%
#chain_arr = np.load(mv.sim_path_MAC + 'test_chain_size/160_mon/chain_6.npy')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

beg=0
end=160

ax.plot(chain_arr[beg-1:beg+1, 0], chain_arr[beg-1:beg+1, 1], chain_arr[beg-1:beg+1, 2],\
        'b--')
ax.plot(chain_arr[beg:end, 0], chain_arr[beg:end, 1], chain_arr[beg:end, 2], 'b-')
ax.plot(chain_arr[end-1:end+1, 0], chain_arr[end-1:end+1, 1], chain_arr[end-1:end+1, 2],\
        'b--')

ax.set_xlabel('x, nm')
ax.set_ylabel('y, nm')
ax.set_zlabel('z, nm')
#plt.show()

xc, yc, zc = np.average(chain_arr, axis=0)

#%% draw sphere
R = 2

u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
x = R * np.cos(u)*np.sin(v) + xc
y = R * np.sin(u)*np.sin(v) + yc
z = R * np.cos(v) + zc
ax.plot_wireframe(x, y, z, color='r')

ax.set_xlim(xc - 3, xc + 3)
ax.set_ylim(yc - 3, yc + 3)
ax.set_zlim(zc - 3, zc + 3)

plt.show()

#%%
def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

angles = []

for i in range(150):
    
    vector_1 = chain_arr[i+1] - chain_arr[i]
    vector_2 = chain_arr[i+2] - chain_arr[i+1]
    
    angles.append(np.rad2deg(angle(vector_1, vector_2)))






