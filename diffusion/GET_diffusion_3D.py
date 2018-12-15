#%% Import
import numpy as np
from matplotlib import pyplot as plt, cm
import os
import importlib
import copy
from itertools import product

import my_functions as mf
import my_variables as mv
import my_indexes as mi
import my_constants as mc
import my_mapping as mm

mf = importlib.reload(mf)
mv = importlib.reload(mv)
mi = importlib.reload(mi)
mc = importlib.reload(mc)
mm = importlib.reload(mm)

os.chdir(mv.sim_path_MAC + 'diffusion')

#%%
mon_rad_matrix = np.load(mv.sim_path_MAC + 'MATRIXES/MATRIX_rad_mon.npy')

#%%
l_xyz = np.array((600, 100, 122))

space = 50
beam_d = 1

x_beg, y_beg, z_beg = (-l_xyz[0]/2, 0, 0)
xyz_beg = np.array((x_beg, y_beg, z_beg))
xyz_end = xyz_beg + l_xyz
x_end, y_end, z_end = xyz_end

step_2nm = 2

x_bins_2nm = np.arange(x_beg, x_end + 1, step_2nm)
y_bins_2nm = np.arange(y_beg, y_end + 1, step_2nm)
z_bins_2nm = np.arange(z_beg, z_end + 1, step_2nm)

bins_2nm = x_bins_2nm, y_bins_2nm, z_bins_2nm

x_grid_2nm = (x_bins_2nm[:-1] + x_bins_2nm[1:]) / 2
y_grid_2nm = (y_bins_2nm[:-1] + y_bins_2nm[1:]) / 2
z_grid_2nm = (z_bins_2nm[:-1] + z_bins_2nm[1:]) / 2

resist_shape = len(x_grid_2nm), len(y_grid_2nm), len(z_grid_2nm)

#%%
D = 2
dx = 2
dy = 2

dt = 0.01

x = x_grid_2nm
z = z_grid_2nm

d_coord = dx

#%%
u = copy.deepcopy(mon_rad_matrix)

add_x = np.zeros((1, len(u[0]), len(u[0, 0])))
u = np.concatenate((add_x, u, add_x), axis=0)

add_y = np.zeros((len(u), 1, len(u[0, 0])))
u = np.concatenate((add_y, u, add_y), axis=1)

add_z = np.zeros((len(u), len(u[0]), 1))
u = np.concatenate((add_z, u, add_z), axis=2)

#%%
nt = 10000

for n in range(nt):
    
    mf.upd_progress_bar(n, nt)
    
    un = copy.deepcopy(u)
    
    u[:, :, -1] = u[:, :, -2]
    
    u[1:-1, 1:-1, 1:-1] = un[1:-1, 1:-1, 1:-1] + D * dt / d_coord**2 * (
        un[2:, 1:-1, 1:-1] - 2 * un[1:-1, 1:-1, 1:-1] + un[0:-2, 1:-1, 1:-1] +
        un[1:-1, 2:, 1:-1] - 2 * un[1:-1, 1:-1, 1:-1] + un[1:-1, 0:-2, 1:-1] +
        un[1:-1, 1:-1, 2:] - 2 * un[1:-1, 1:-1, 1:-1] + un[1:-1, 1:-1, 0:-2]
        )

#%%
fig = plt.figure()
ax = fig.gca(projection='3d')

X, Z = np.meshgrid(x, z)

surf = ax.plot_surface(X, Z, np.sum(u[1:-1, 1:-1, 1:-1], axis=1).transpose(), rstride=1,\
    cstride=1, cmap=cm.viridis, linewidth=0, antialiased=True)
#ax.set_zlim(0, 2.5)
ax.set_xlabel('$x$')
ax.set_ylabel('$z$')
