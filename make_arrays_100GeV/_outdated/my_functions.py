from random import random
import numpy as np
import sys
from scipy import interpolate
import my_arrays as ma
import my_variables as mv
import matplotlib.pyplot as plt
import my_array_joiner as ma

#%% Non-simulation functions
def upd_progress_bar(total, progress):
    barLength, status = 20, ''
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, '\r\n'
    block = int(round(barLength * progress))
    text = '\r[{}] {:.0f}% {}'.format(
        '#' * block + '-' * (barLength - block), round(progress * 100, 0),\
        status)
    sys.stdout.write(text)
    sys.stdout.flush()

def log_interp1d(xx, yy, kind='linear'):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interpolate.interp1d(logx, logy, kind=kind)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

#%% Beam functions
def get_norm_density(x, mu, sigma):
    y = 1/(sigma*np.sqrt(2*np.pi))*np.exp((-1)*(x - mu)**2/(2*sigma**2))
    return y

def get_x_y_round_beam(x0, y0, R):
    rho = random()*R
    phi = 2*np.pi*random()
    x = x0 + rho*np.cos(phi)
    y = y0 + rho*np.sin(phi)
    return x, y

def get_x_y_square_beam(x0, y0, D):
    x = x0 + D*random()
    y = y0 + D*random()
    return x, y

#%% Simulation functions
def get_closest_el_ind(array, val):
    ind = np.argmin(np.abs(array - val))
    return ind
   
def get_closest_norm_int_ind(array, value):
    for i in range(len(array) - 1):           
        if array[i] <= value <= array[i + 1]:
            return i
            break
    return len(array)

def get_atom_id(d_PMMA, E_ind, z):
    if d_PMMA == 0:
        return 3
    else:
        if 0 <= z <= d_PMMA: # H*8 + C*5 + O*2
            array = ma.PMMA_ATOM_CS_SUM_NORM[int(E_ind), :]
            atom_id = get_closest_norm_int_ind(array, random())
            return atom_id
        elif z > d_PMMA:
            return 3
        else:
            return 4

def get_mean_free_path(atom_id, E_ind):
    mfp = 1/(mv.CONC[atom_id]*np.sum(ma.CS_TOTAL[atom_id][E_ind, :]))
    return mfp
    
def get_collision_id(atom_id, E_ind):
    array = ma.CS_TOTAL_SUM_NORM[atom_id][E_ind, :]
    coll_id = get_closest_norm_int_ind(array, random())
    return coll_id

def get_On(atom_id, coll_id, E_ind, O_prev):
    if coll_id == 0:
        theta_arr = ma.DIFF_CS_INT_NORM[atom_id][E_ind, :]
        theta_ind = get_closest_norm_int_ind(theta_arr, random())
        On = get_O_matrix(2*np.pi*random(), ma.theta_arr[theta_ind], O_prev)
        return On
    else:
        return O_prev

def get_O_matrix(phi, theta, O_prev):
    Wn = np.mat([[      np.cos(phi),                np.sin(phi),             0],\
        [-np.sin(phi)*np.cos(theta),  np.cos(phi)*np.cos(theta), np.sin(theta)],\
        [ np.sin(phi)*np.sin(theta), -np.cos(phi)*np.sin(theta), np.cos(theta)]])
    return Wn*O_prev

def get_dxdydz(d_PMMA, mfp, On, z):
    x0 = np.mat([[0], [0], [1]])
    ds = -mfp * np.log(random())
    dxdydz = On.transpose()*x0*ds
    if d_PMMA > 0 and (z - d_PMMA) * (z + dxdydz[2] - d_PMMA) < 0:
        dxdydz *= (d_PMMA - z)/np.linalg.norm(dxdydz)
    return dxdydz

def get_final_On_and_O2nd(E_prev, E2nd, On):
    phi_ion = 2*np.pi*random()
    omega, t = E2nd/E_prev, E_prev/mv.m_eV
    alpha = np.arcsin(np.sqrt(2*omega/(2 + t - t*omega)))
    gamma = np.arcsin(np.sqrt(2*(1 - omega)/(2 + t*omega)))
    O_final = get_O_matrix(phi_ion, alpha, On)
    O2nd = get_O_matrix(np.pi + phi_ion, gamma, On)
    return O_final, O2nd

def get_coll_data(d_PMMA, E_prev, O_prev, x, z):
    is2nd = False
    E2nd = 0
    O2nd = O_prev*0
    
    E_ind = get_closest_el_ind(ma.E_arr, E_prev)
    atom_id = get_atom_id(d_PMMA, z, E_ind)
    mfp = get_mean_free_path(atom_id, E_ind)
    coll_id = get_collision_id(atom_id, E_ind)
    On = get_On(atom_id, coll_id, E_ind, O_prev)
    dxdydz = get_dxdydz(d_PMMA, mfp, On, z)
    
    if z == 0 or np.abs(z - d_PMMA) < 1e-10:
        return (atom_id, coll_id, E_prev, dxdydz.transpose(), O_prev, False, 0, 0, 0)
    
#    if coll_id < 2:
#        dE = ma.DE_FULL[atom_id][E_ind, coll_id]
#        E = E_prev - dE
    if coll_id == 0:
        dE = 0
        E = E_prev
    
    elif coll_id == 1:
        dE = ma.DE_EXC[atom_id][E_ind]
        E = E_prev - dE
        
    else: # IONIZATION
        subshell_id = coll_id - 2
        
        spectra_line = ma.ION_SPECTRA[atom_id][subshell_id][E_ind, :]
        E_ext = ma.ION_E_EXT[atom_id][subshell_id]
        E2nd = E_ext[get_closest_el_ind(spectra_line, random())]
        
        dE = 0
#        dE = E2nd + ma.ION_E_BIND[atom_id][subshell_id]
        
        if 15 <= E2nd <= E_prev:
            is2nd = True
            On, O2nd = get_final_On_and_O2nd(E_prev, E2nd, On)
#        dE = E2nd + ma.ION_E_BIND[atom_id][E_ind, subshell_id]
            dE = E2nd + ma.ION_E_BIND[atom_id][subshell_id]
            
        E = E_prev - dE
        
#    if dE > 2e+4:
#        print(atom_id, coll_id, E, dxdydz.transpose(), On, is2nd, E2nd, O2nd, dE)
        
    return atom_id, coll_id, E, dxdydz.transpose(), On, is2nd, E2nd, O2nd, dE

def get_TASKS_and_sim_data(d_PMMA, TASKS, tr_num, par_num, E0, x0y0z0, O_prev):
    # DATA array structure:
    # track_num | parent_num | aid | pid | E | x | y | z | dE
    E = E0
#    sim_data = np.zeros((int(2e+4), 9))*np.nan
    sim_data = np.zeros((int(2e+7), 9))*np.nan
    pos = 0
    sim_data[pos, :] = np.hstack((tr_num, par_num, np.nan, np.nan, E0, x0y0z0, np.nan))
    
    while E > 15:
        x = sim_data[pos, 5]
        z = sim_data[pos, 7]
        if z < 0:
            break
        aid, cid, E, dxdydz, On, is2nd, E2nd, O2nd, dE =\
            get_coll_data(d_PMMA, E, O_prev, x, z)
        new_task = []
        if (is2nd):
            new_task = [tr_num, E2nd, sim_data[pos, 5:-1], O2nd]
            TASKS.append(new_task)
        sim_data[pos, 2] = aid
        sim_data[pos, 3] = cid
        sim_data[pos, 8] = dE
        sim_data[pos + 1, :] = np.concatenate(([[tr_num]], [[par_num]], [[np.nan]],\
                 [[np.nan]], [[E]], sim_data[pos, 5:-1] + dxdydz, [[np.nan]]), axis=1)
        O_prev = On
        pos += 1
    
    sim_data = np.delete(sim_data, np.where(np.isnan(sim_data[:, 0])), axis=0)    
    return TASKS, sim_data

def create_TASKS(E0, n_tracks):
    O0 = np.eye(3)
    TASKS = [None]*n_tracks
    for i in range(n_tracks):
        x0, y0 = 0, 0
        coords = np.array(np.hstack((x0, y0, 0)))
        task = [np.nan, E0, coords, O0]
        TASKS[i] = task
    return TASKS

def get_DATA(beam_pos, D, d_PMMA, n_tracks):
#    n_coords = int(5e+3)
    n_coords = int(1e+3)
#    E0 = 20e+3
    E0 = 1e+6
    TASKS = create_TASKS(E0, n_tracks)
    DATA = np.zeros((n_coords*n_tracks, 9))*np.nan
    dataline_pos = 0
    track_num = 0
    
    # create DATA file for TASKS
    while track_num < len(TASKS):
        upd_progress_bar(len(TASKS), track_num + 1)
        task = TASKS[track_num]
        par_num, E0, coords, O0 = task[0], task[1], task[2], task[3]
        TASKS, tr_data = get_TASKS_and_sim_data(d_PMMA, TASKS, track_num,\
                                                par_num, E0, coords, O0)
        DATA[dataline_pos:dataline_pos + len(tr_data), :] = tr_data
        dataline_pos += len(tr_data)
        track_num += 1

    DATA = np.delete(DATA, np.where(np.isnan(DATA[:, 2])), axis=0)
    
    return DATA

#%% Plot DATA
def plot_DATA(DATA):
    fig, ax = plt.subplots()
    for tn in range(int(np.max(DATA[:, 0]))):
        if len(np.where(DATA[:, 0] == tn)[0]) == 0:
            continue
        beg = np.where(DATA[:, 0] == tn)[0][0]
        end = np.where(DATA[:, 0] == tn)[0][-1] + 1
        ax.plot(DATA[beg:end, 5], DATA[beg:end, 7], linewidth=0.7)
    #    inds_el = beg + np.where(DATA[beg:end, 3] == 0)[0]
    #    inds_ion = beg + np.where(DATA[beg:end, 3] >= 2)[0]
    #    inds_exc = beg + np.where(DATA[beg:end, 3] == 1)[0]
    #    ax.plot(DATA[inds_el, 5], DATA[inds_el, 7], 'r.')
    #    ax.plot(DATA[inds_ion, 5], DATA[inds_ion, 7], 'b.')
    #    ax.plot(DATA[inds_exc, 5], DATA[inds_exc, 7], 'g.')
    
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title('Direct Monte-Carlo simulation')
    plt.xlabel('x, cm')
    plt.ylabel('z, cm')
    plt.axis('on')
    plt.grid('on')
    plt.gca().invert_yaxis()
    plt.show()