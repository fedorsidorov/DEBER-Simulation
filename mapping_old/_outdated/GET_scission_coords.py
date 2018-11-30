import os
import numpy as np
import _my_functions as mf

#%%
eps = 1e-3

n_cubes_in_layer = 100
cells_in_cube = 20, 20, 20

cell_size = 0.5
#monomer_size = 0.25

step = cell_size

chain_dir = '9_jul' + os.sep + 'HOME_chains_sorted'
e_data_dir = '9_jul' + os.sep + 'cubes_1C_L3_ion'
dest_dir = '9_jul' + os.sep + 'results_1C_L3_ion'

if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

#%%
chain_filenames = os.listdir(chain_dir)

n_file = 0

for chain_fname in chain_filenames:
    
    if chain_fname == '.DS_Store':
        continue
    
    mf.upd_progress_bar(n_file, len(chain_filenames))
    
    ## load chain array
    ch_raw_arr = np.load(chain_dirname + os.sep + chain_fname)
    chain_L = len(ch_raw_arr)
    
    ## make chain shift to get positive coordinates
    ch_arr = ch_raw_arr - ch_raw_arr.min(axis=0)
    
    ## get number of cells required
    nc_x, nc_y, nc_z = (ch_arr[:, :-1].max(axis=0) // cells_in_cube + 1).astype(int)
    
    if (nc_x > 2) or (nc_y > 2) or (nc_z > 4):
        print('Oversize')
        n_file += 1
        continue
    
    ## array for future e data
    e_data = np.zeros((0, 3))
    
    ## set z shift to chain
    z_shift = mf.choice(list(range(4)))
    
    ## get array with electron events
    for k in range(nc_z):                
        
        ## get x-y cubes numbers for every z layer
        cubes_xy_nums = mf.choice(np.arange(n_cubes_in_layer), nc_x*nc_y, replace=False)
        xy_pos = 0
        
        for j in range(nc_y):
            for i in range(nc_x):
            
                z_pos = z_shift + k
                
                if z_pos > 3:
                    z_pos = 3
                
                ## load required cube
                cube_events_raw = np.load(e_data_dirname + os.sep + 'cube_' +\
                    str(z_pos) + '_' + str(cubes_xy_nums[xy_pos]) + '.npy')
                
                ## convert cube indexes to nanometers
                cube_events = (cube_events_raw + np.array((i, j, k)) * cells_in_cube)\
                    * cube_cell_size
                
                ## add cube events to e_data
                e_data = np.vstack((e_data, cube_events))
                
                ## increment xy_pos
                xy_pos += 1
    
    result_coll_arr = np.zeros((len(e_data), 4))*np.nan
    coll_pos = 0

    beg = 0
    set_beg = False

    for e_line in e_data:

        x_e, y_e, z_e = e_line
        
        ## find collisions with e-data in chain array
        statements = (ch_arr[:, 0] >= x_e - eps,\
                      ch_arr[:, 0] <= x_e + cube_cell_size + eps,\
                      ch_arr[:, 1] >= y_e - eps,\
                      ch_arr[:, 1] <= y_e + cube_cell_size + eps,\
                      ch_arr[:, 2] >= z_e - eps,\
                      ch_arr[:, 2] <= z_e + cube_cell_size + eps)
    
        ch_shift_coll_arr = ch_arr[np.where(np.logical_and.reduce(statements))]

        ## if there is 1 or more collisions
        if len(ch_shift_coll_arr) > 0:
            ind = mf.choice(tuple(range(len(ch_shift_coll_arr))))
            result_coll_arr[coll_pos, :] = ch_shift_coll_arr[ind, :]
            coll_pos += 1

    result = mf.delete_nan_rows(result_coll_arr)
    fname = 'Ch_' + str(n_file) + '_L_' + str(chain_L) + '.npy'
    np.save(dest_dirname + os.sep + fname, result)
    print(fname)
    
    n_file += 1
