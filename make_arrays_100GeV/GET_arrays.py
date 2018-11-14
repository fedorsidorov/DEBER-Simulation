import numpy as np
import sys

#import matplotlib.pyplot as plt
#import my_arrays as ma
#import my_variables as mv

sys.path.append('../MODULES')
import my_functions as mf

#E_arr = np.power(10, np.arange(1, 4.3, 0.01))
E_arr_x = np.power(10, np.arange(1, 11, 0.01))
theta_arr = np.hstack((np.power(10, np.arange(-2, np.log10(np.pi), 0.01)), np.pi))

#np.save('IND_ARRAYS/E_arr.npy', E_arr)
#np.save('IND_ARRAYS/E_arr_x.npy', E_arr_x)
#np.save('IND_ARRAYS/theta_arr.npy', theta_arr)

#%% Elastic

def get_elastic_total_cs(el_name):
    
    path = el_name.join(['', '_data/Elastic/', '_10_0.dat'])
    raw_arr = np.loadtxt(path, delimiter=' ', usecols=(1, 2))
    
    inds = np.where(E_arr_x < raw_arr[0, 0]*1e+6)[0]
    add_arr = np.hstack((E_arr_x[inds].reshape((len(inds), 1))/1e+6, np.zeros((len(inds), 1))))
    raw_arr = np.vstack((add_arr, raw_arr))
    
    result = mf.log_interp1d(raw_arr[:, 0]*1e+6, raw_arr[:, 1]*1e-24)(E_arr_x)
    result[np.where(np.isnan(result))] = 0

    return result


def get_elastic_diff_cs(el_name):

    path = el_name.join(['', '_data/Elastic/', '_8_22.dat'])
    raw_arr = np.loadtxt(path, delimiter=' ', usecols=(1, 2, 3))
    
    arr = np.zeros(np.shape(raw_arr))*np.NaN
    arr[:, 0] = raw_arr[:, 0]*1e+6
    arr[:, 1] = np.arccos(-raw_arr[:, 1] + 1)
    arr[:, 2] = raw_arr[:, 2]
    
    E_unq = np.unique(arr[:, 0])
    
    # интерполяция по углу
    diff_arr_0 = np.zeros((len(E_unq), len(theta_arr)))*np.NaN
    
    for i in range(len(E_unq)):
        inds = np.where(arr[:, 0] == E_unq[i])
        diff_arr_0[i, :] = mf.log_interp1d(arr[inds, 1][0, :], arr[inds, 2][0, :])(theta_arr)
    
    # интерполяция по энергии
    diff_arr_1 = np.zeros((len(E_arr_x), len(theta_arr)))*np.NaN
    
    for j in range(len(theta_arr)):
        diff_arr_1[:, j] = mf.log_interp1d(E_unq, diff_arr_0[:, j])(E_arr_x)
    
    # получим полные сечения
    elastic_total_cs = get_elastic_total_cs(el_name)
    
    # зададим правильную нормировку
    diff_arr_2 = np.zeros(np.shape(diff_arr_1))
    
    for i in range(len(E_arr_x)):
        integral = np.trapz(2*np.pi*diff_arr_1[i, :]*np.sin(theta_arr), x=theta_arr)
        diff_arr_2[i, :] = diff_arr_1[i, :] / integral * elastic_total_cs[i]
    
    # и наконец заполним массивы с дифференциальными сечениями
    diff_arr_int = np.zeros((len(E_arr_x), len(theta_arr)))
    diff_arr_int_norm = np.zeros((len(E_arr_x), len(theta_arr)))
    
    for i in range(len(E_arr_x)):
        s_tot = np.trapz(2*np.pi*diff_arr_2[i,:]*np.sin(theta_arr), x=theta_arr)
    
        for j in range(len(theta_arr) - 1):
            diff_arr_int[i, j] = np.trapz(2*np.pi*diff_arr_2[i, 0:j + 1]*\
                       np.sin(theta_arr[0:j + 1]), x=theta_arr[0:j + 1])
            diff_arr_int_norm[i, j] = diff_arr_int[i, j] / s_tot
        
        diff_arr_int[i, len(theta_arr) - 1] = s_tot
        diff_arr_int_norm[i, len(theta_arr) - 1] = 1

    return diff_arr_int_norm

a = get_elastic_diff_cs('Si')


#%% Excitation

def get_excitation_total_cs(el_name):
    
    path = el_name.join(['', '_data/Excitation/', '_83_0.dat'])
    raw_arr = np.loadtxt(path, delimiter=' ', usecols=(1, 2))
    
    inds = np.where(E_arr_x < raw_arr[0, 0]*1e+6)[0]
    add_arr = np.hstack((E_arr_x[inds].reshape((len(inds), 1))/1e+6, np.zeros((len(inds), 1))))
    raw_arr = np.vstack((add_arr, raw_arr))
    
    result = mf.log_interp1d(raw_arr[:, 0]*1e+6, raw_arr[:, 1]*1e-24)(E_arr_x)
    result[np.where(np.isnan(result))] = 0

    return result


def get_excitation_dE(el_name):
    
    path = el_name.join(['', '_data/Excitation/', '_83_11.dat'])
    raw_arr = np.loadtxt(path, delimiter=' ', usecols=(1, 2))
    
    inds = np.where(E_arr_x < raw_arr[0, 0]*1e+6)[0]
    add_arr = np.hstack((E_arr_x[inds].reshape((len(inds), 1))/1e+6, np.zeros((len(inds), 1))))
    raw_arr = np.vstack((add_arr, raw_arr))
    
    result = mf.log_interp1d(raw_arr[:, 0]*1e+6, raw_arr[:, 1]*1e+6)(E_arr_x)
    result[np.where(np.isnan(result))] = 0
    
    return result


#%% Ionization

def get_ionization_cs(el_name, ss_name):
    
    path_0 = el_name.join(['', '_data/Ionization/'])
    path_1 = ss_name.join(['', '/0_81_0_91.dat'])
    raw_arr = np.loadtxt(path_0 + path_1, delimiter=' ', usecols=(1, 2))
    
    inds = np.where(E_arr_x < raw_arr[0, 0]*1e+6)[0]
    add_arr = np.hstack((E_arr_x[inds].reshape((len(inds), 1))/1e+6, np.zeros((len(inds), 1))))
    raw_arr = np.vstack((add_arr, raw_arr))
    
    result = mf.log_interp1d(raw_arr[:, 0]*1e+6, raw_arr[:, 1]*1e-24)(E_arr_x)
    result[np.where(np.isnan(result))] = 0
    
    return result


def get_ionization_E2nd(el_name, ss_name):
    
    path_0 = el_name.join(['', '_data/Ionization/'])
    path_1 = ss_name.join(['', '/19_81_10_91.dat'])
    raw_arr = np.loadtxt(path_0 + path_1, delimiter=' ', usecols=(1, 2))
    
    inds = np.where(E_arr_x < raw_arr[0, 0]*1e+6)[0]
    add_arr = np.hstack((E_arr_x[inds].reshape((len(inds), 1))/1e+6, np.zeros((len(inds), 1))))
    raw_arr = np.vstack((add_arr, raw_arr))
    
    result = mf.log_interp1d(raw_arr[:, 0]*1e+6, raw_arr[:, 1]*1e+6)(E_arr_x)
    result[np.where(np.isnan(result))] = 0

    return result


def get_ionization_2nd_spectra(el_name, ss_name, ion_E2nd):

    path_0 = el_name.join(['', '_data/Ionization/'])
    path_1 = ss_name.join(['', '/19_81_21_91.dat'])
    file_raw = np.loadtxt(path_0 + path_1, delimiter=' ', usecols=(1, 2, 3))
    
    # конветируем все в эВ
    file_0 = np.concatenate((file_raw[:, 0][:, None]*1e+6, file_raw[:, 1][:, None]*1e+6,\
                      file_raw[:, 2][:, None]*1e-6), axis=1)
    E_min, E_max = np.min(file_0[:, 1]), np.max(file_0[:, 1])
    
    E_ext_pre = np.power(10, np.arange(np.log10(E_min) + 0.01, np.log10(E_max) - 0.01, 0.01))
    E_ext = np.hstack((E_min, E_ext_pre, E_max))
    Es = len(E_arr_x)
    Exs = len(E_ext)
    
    
    E_unq = np.unique(file_0[:, 0])
    
    E_end = E_unq[-1]
    block_end = file_0[np.where(file_0[:, 0] == E_end)]
    beg = np.where(block_end[:, 1] > 1e+4)[0][0] # may vary
    end = np.where(block_end[:, 1] > 1e+8)[0][0] # may vary
    x = block_end[beg:end, 1]
    y = block_end[beg:end, 2]
    A = np.vstack([np.log10(x), np.ones(len(x))]).T
    k, b = np.linalg.lstsq(A, np.log10(y))[0]
    
    E_ext_y = 10**b * E_ext**k
    
    # Расширим данные
    file_1 = np.zeros((0, 3))
    for e in E_unq:
        block = file_0[np.where(file_0[:, 0] == e), :][0]
        if block[0, 1] > E_min:
            block = np.vstack((np.array((e, E_min, block[0, 2])), block))
        inds_x = np.where(E_ext >= block[-1, 1])[0]
        inds_y = np.where(E_ext_y < block[-1, 2])[0]
        if inds_x[0] > inds_y[0]:
            inds = inds_x
        else:
            inds = inds_y
        add_arr = np.concatenate((np.ones(np.shape(E_ext[inds]))[:, None]*e, E_ext[inds][:, None],\
                                 E_ext_y[inds][:, None]), axis=1)
        if np.shape(block)[0] >= 10: # and e != E_unq[-1]: # to prevent sharp transitions
            file_1 = np.vstack((file_1, block[:-6, :], add_arr))
        else:
            file_1 = np.vstack((file_1, block, add_arr))
    
    # Добавим наименьшую энергию, если нужно
    E_unq_x = np.copy(E_unq)
    file_1_x = np.copy(file_1)
    
    if E_unq_x[0] > E_arr_x[0]:
        E_unq_x = np.hstack((E_arr_x[0], E_unq))
        pre_block = file_1[np.where(file_1[:, 0] == E_unq_x[1])]
        pre_block[:, 0] = E_arr_x[0]
        file_1_x = np.vstack((pre_block, file_1))

    spectra_0 = np.zeros((np.shape(E_unq_x)[0], Exs))*np.nan
    
    for i in range(np.shape(E_unq_x)[0]):
        inds = np.where(file_1_x[:, 0] == E_unq_x[i])
        spectra_0[i, :] = mf.log_interp1d(file_1_x[inds, 1][0, :],\
                   file_1_x[inds, 2][0, :])(E_ext)
    
    #
    spectra_1 = np.zeros((Es, Exs))*np.nan
    
    for j in range(np.shape(E_ext)[0]):
        spectra_1[:, j] = mf.log_interp1d(E_unq_x, spectra_0[:, j])(E_arr_x)
    
    # Cut data - mean should be equal to ion_E2nd_files[i]
    spectra_2 = np.zeros(np.shape(spectra_1))
    
    for i in range(Es):
        for j in range(Exs):
            if np.trapz(spectra_2[i, :j]*E_ext[:j], x=E_ext[:j]) <= ion_E2nd[i]:
                spectra_2[i, j] = spectra_1[i, j]
            else:
                break
    
    ion_spectra_int = np.zeros((Es, Exs))
    ion_spectra_int_norm = np.zeros((Es, Exs))
    
    for i in range(Es):
        sum_tot = np.trapz(spectra_2[i, :], x=E_ext)
        
        for j in range(Exs - 1):
            ion_spectra_int[i, j] = np.trapz(spectra_2[i, 0:j + 1], x=E_ext[0:j + 1])
            ion_spectra_int_norm[i, j] = ion_spectra_int[i, j]/sum_tot
            
        ion_spectra_int[i, Exs - 1] = sum_tot
        ion_spectra_int_norm[i, Exs - 1] = 1
    
    return ion_spectra_int_norm


#el_name, ss_name = 'Si', 'K'
#ion_E2nd = get_ionization_E2nd(el_name, ss_name)
#a = get_ionization_2nd_spectra(el_name, ss_name, ion_E2nd)
