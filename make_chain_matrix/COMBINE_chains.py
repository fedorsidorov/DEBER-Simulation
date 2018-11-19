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
source_dir = mv.sim_path_MAC + 'CHAINS/CHAINS_950K_122nm_10k/'
dest_dir = mv.sim_path_MAC + 'CHAINS/CHAINS_950K_122nm_10k_shifted/'

N_chains = 10000

chain_bank = [[] for i in range(N_chains)]
L_arr = np.zeros(N_chains)

## load chains into bank
for i in range(N_chains):
    
    mf.upd_progress_bar(i, N_chains)
    
    chain_bank[i] = np.load(source_dir + 'chain_' + str(i % len(chain_bank)) + '.npy')
    L_arr[i] = len(chain_bank[i])
    
#%% check log_mw
log_mw = np.log10(L_arr * 100)
plt.hist(log_mw)
plt.title('Chain mass distribution')
plt.xlabel('log(m$_w$)')
plt.ylabel('probability')
plt.legend()
plt.grid()
plt.show()

#%% prepare histograms
l_xyz = np.array((100, 100, 122))
space = 100

x_beg, y_beg, z_beg = (0, 0, 0)
xyz_beg = np.array((x_beg, y_beg, z_beg))
xyz_end = xyz_beg + l_xyz
x_end, y_end, z_end = xyz_end

step_prec = 10
step_mono = 0.28

bins_total = np.array(np.hstack((xyz_beg.reshape(3, 1), xyz_end.reshape(3, 1))))

#%%
x_bins_prec = np.arange(x_beg, x_end + 1, step_prec)
y_bins_prec = np.arange(y_beg, y_end + 1, step_prec)
z_bins_prec = np.arange(z_beg, z_end + 1, step_prec)

bins_prec = [x_bins_prec, y_bins_prec, z_bins_prec]

hist_total = np.zeros((1, 1, 1))
hist_prec = np.zeros((len(x_bins_prec) - 1, len(y_bins_prec) - 1, len(z_bins_prec) - 1))

# create chain_list and check density
chain_list = []

V = np.prod(l_xyz) * (1e-7)**3 ## cm^3
m_mon = 1.66e-22 ## g
rho = 1.19 ## g / cm^3

i = 0

while True:
    
    if i % 1000 == 0:
        print(i, 'chains are added')
    
    if np.sum(hist_total) * m_mon / V >= rho:
        print('Needed density is achieved')
        break
    
    else:
        
        ii = i % len(chain_bank)
        new_chain = chain_bank[ii]
        
        x_shift = mf.uniform(-space, l_xyz[0] + space)
        y_shift = mf.uniform(-space, l_xyz[0] + space)
        
        new_chain_shift = new_chain + np.array((x_shift, y_shift, 0))
        
        if new_chain_shift.max(axis=0)[0] < x_beg or new_chain_shift.max(axis=0)[1] < y_beg\
        or new_chain_shift.min(axis=0)[0] > x_end or new_chain_shift.min(axis=0)[1] > y_end:
            continue
        
        chain_list.append(new_chain_shift)
        
        hist_total += np.histogramdd(new_chain_shift, bins=bins_total)[0]
        hist_prec += np.histogramdd(new_chain_shift, bins=bins_prec)[0]
        
    i += 1

#%% check density
density_total = hist_total[0][0][0] * m_mon / V
density_prec = hist_prec * m_mon / V * (len(x_bins_prec) - 1) * (len(y_bins_prec) - 1) *\
    (len(z_bins_prec) - 1)

#%% save chains to files
i = 0

for chain in chain_list:
    
    mf.upd_progress_bar(i, len(chain_list))
    np.save(dest_dir + 'chain_shift_' + str(i) + '.npy', chain)
    i += 1

#%% check n_mon in 0.28 nm cube
#array_mono = np.reshape(hist_mono, (np.prod(np.shape(hist_mono)),))
#array_mono_hist = np.histogram(array_mono, bins=10)[0]
#mf.print_histogram(array_mono, user_bins=10, is_normed=False)

#%% cut chains to cube shape
chain_cut_list = []

for chain in chain_list:
    
    statements = [chain[:, 0] >= x_beg, chain[:, 0] <= x_end,
                  chain[:, 1] >= y_beg, chain[:, 1] <= y_end]
    inds = np.where(np.logical_and.reduce(statements))[0]
    
    beg = 0
    end = -1
    
    for i in range(len(inds) - 1):
        if inds[i+1] > inds[i] + 1 or i == len(inds) - 2:
            end = i + 1
            chain_cut_list.append(chain[inds[beg:end], :])
            beg = i + 1

#%% get nice 3D picture
l_xy = len(np.arange(0, 101, 1))
l_z = len(np.arange(0, 123, 1))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for chain in chain_list[0:-1:50]:
    ax.plot(chain[:, 0], chain[:, 1], chain[:, 2])

ax.plot(np.arange(0, 101, 1), np.ones(l_xy)*0, np.ones(l_xy)*0, 'k')
ax.plot(np.arange(0, 101, 1), np.ones(l_xy)*100, np.ones(l_xy)*0, 'k')
ax.plot(np.arange(0, 101, 1), np.ones(l_xy)*0, np.ones(l_xy)*122, 'k')
ax.plot(np.arange(0, 101, 1), np.ones(l_xy)*100, np.ones(l_xy)*122, 'k')

ax.plot(np.ones(l_xy)*0, np.arange(0, 101, 1), np.ones(l_xy)*0, 'k')
ax.plot(np.ones(l_xy)*100, np.arange(0, 101, 1), np.ones(l_xy)*0, 'k')
ax.plot(np.ones(l_xy)*0, np.arange(0, 101, 1), np.ones(l_xy)*122, 'k')
ax.plot(np.ones(l_xy)*100, np.arange(0, 101, 1), np.ones(l_xy)*122, 'k')

ax.plot(np.ones(l_z)*0, np.ones(l_z)*0, np.arange(0, 123, 1), 'k')
ax.plot(np.ones(l_z)*100, np.ones(l_z)*0, np.arange(0, 123, 1), 'k')
ax.plot(np.ones(l_z)*0, np.ones(l_z)*100, np.arange(0, 123, 1), 'k')
ax.plot(np.ones(l_z)*100, np.ones(l_z)*100, np.arange(0, 123, 1), 'k')

plt.xlim(x_beg, x_end)
plt.ylim(y_beg, y_end)
plt.title('Polymer chain simulation')
ax.set_xlabel('x, nm')
ax.set_ylabel('y, nm')
ax.set_zlabel('z, nm')
plt.show()