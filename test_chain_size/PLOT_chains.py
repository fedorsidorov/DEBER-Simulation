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
chain_arr = np.load(mv.sim_path_MAC + 'test_chain_size/160_mon/chain_6.npy')

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

## draw sphere
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
