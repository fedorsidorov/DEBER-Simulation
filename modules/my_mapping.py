#%% Import
import numpy as np
import importlib
import my_functions as mf
import my_variables as mv

mf = importlib.reload(mf)
mv = importlib.reload(mv)

#%% resist_matrix
n_chain_ind = 0
beg_mon, mid_mon, end_mon = 0, 1, 2
free_mon, free_rad_mon = 10, 20

#%% chain_table
x_ind, y_ind, z_ind = 0, 1, 2
mon_line_pos_ind = 3
mon_type_ind = -1
uint16_max = 65535

#%% constants
N_chains_total = 63306
N_mon_chain_max = 9780
N_mon_cell_max = 810

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
## changes monomer type
def rewrite_mon_type(resist_matrix, chain_table, n_chain, n_mon, new_type):
    ## change mon_type in chain_inv_matrix
    chain_table[n_chain, n_mon, mon_type_ind] = new_type
    ## define x, y, z of monomer
    x_ind, y_ind, z_ind, mon_line_pos = chain_table[n_chain, n_mon, :mon_type_ind]
    ## if we aren't outside the resist area of interest
    if not x_ind == y_ind == z_ind == mon_line_pos == uint16_max:
        resist_matrix[x_ind, y_ind, z_ind, mon_line_pos, mon_type_ind] = new_type


## choose one of existing particles to interact electron with
def get_resist_part_ind(resist_matrix, x_ind, y_ind, z_ind):
    ## indexes of existing particle lines
    resist_part_inds = np.where(resist_matrix[x_ind, y_ind, z_ind, :, n_chain_ind] !=\
                         uint16_max)[0]
    ## if no free particles
    if len(resist_part_inds) == 0:
            return -1
    return mf.choice(resist_part_inds)


## choose of monomer type
def get_mon_type():
    return mf.choice([0, 2])


## get neede chain_table line
def get_chain_table_line(chain_table, n_chain, n_mon):
    return chain_table[n_chain, n_mon]


## get n events in the cell
def get_n_events(e_matrix, x_ind, y_ind, z_ind):
    return e_matrix[x_ind, y_ind, z_ind]


## get particle line from 
def get_resist_part_line(resist_matrix, x_ind, y_ind, z_ind, resist_part_ind):
    return resist_matrix[x_ind, y_ind, z_ind, resist_part_ind, :]


## convert monomer type to monomer kind
def mon_type_to_kind(mon_type):
    ## with ester group
    if mon_type in [-1, 0, 1]:
        return mon_type
    ## W/O ester group
    else:
        return mon_type - 10


## convert 65536 to -1
def correct_mon_type(mon_type):
    if mon_type == uint16_max:
        return -1
    return mon_type


## calculate local AVG chain length distribution
def get_local_chain_len(chain_table):
    
    chain_sum_len_matrix = np.zeros(resist_shape)
    n_chains_matrix = np.zeros(resist_shape)
    
    for idx, chain in enumerate(chain_table):
        
        mf.upd_progress_bar(idx, N_chains_total)
        
        beg_ind = 0
        
        while True:
            
            if beg_ind >= N_mon_chain_max:
                break
                        
            if chain[beg_ind, mon_type_ind] in [2, 12]:
                beg_ind += 1
                continue
            
            if chain[beg_ind, mon_type_ind] not in [uint16_max, 9]:
                print('mon_type', chain[beg_ind, mon_type_ind])
                print('chain indexing error!')
            
            st_1 = chain[beg_ind:, mon_type_ind] == 1
            st_2 = chain[beg_ind:, mon_type_ind] == 11
            
            where_result = np.where(np.logical_or(st_1, st_2))[0]
            
            if len(where_result) == 0:
                break
            
            end_ind = beg_ind + where_result[0]
            now_chain_len = end_ind - beg_ind
            
            inds_list = []
            
            for mon_line in chain[beg_ind:end_ind+1]:
                
                x_ind, y_ind, z_ind = mon_line[:3]
                
                if x_ind == y_ind == z_ind == uint16_max:
                    continue
                
                now_inds = [x_ind, y_ind, z_ind]
                
                if now_inds in inds_list:
                    continue
                
                chain_sum_len_matrix[x_ind, y_ind, z_ind] += now_chain_len
                n_chains_matrix[x_ind, y_ind, z_ind] += 1
                
                inds_list.append(now_inds)
            
            beg_ind = end_ind + 1
    
    return chain_sum_len_matrix, n_chains_matrix

#%%








