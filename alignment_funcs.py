
import numpy as np 
from particle_utils import string_to_array
import os
import pickle

def theta_dot(theta1, theta2):
    ans = np.sin(theta1)* np.sin(theta2) + np.cos(theta1)*np.cos(theta2)
    return ans
     

def upscore4(thetas):
    score = np.zeros(np.shape(thetas)[1])
    for i in range(6):
        new_score = theta_dot(thetas[i,:], thetas[(i + 1)%6,:])
        score += new_score
 

def hex_align(thetas_dict):
    keys = thetas_dict.keys()
    keys_as_arrays = [string_to_array(key) for key in keys] 
    keys_no_ending = np.unique([arr_key[:-1] for arr_key in keys_as_arrays], axis = 0)
    total_alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]]))
    for block in keys_no_ending:
        hex_keys = np.array([thetas_dict[str(list(np.append(block, [ending])))] for ending in range(6)])
        block_alignment = upscore4(hex_keys)
        total_alignment += block_alignment
    return total_alignment / (len(keys_no_ending) * 6)

def neighbor_align(thetas_dict, neighbor_chart):
    keys = thetas_dict.keys()
    total_alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
    for key in keys:
        alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
        for neighbor in neighbor_chart[key]:
            alignment += theta_dot(thetas_dict[key], thetas_dict[neighbor])
        total_alignment += alignment
    return total_alignment / (len(list(keys)) *2 * len(neighbor_chart[list(keys)[0]])) # * 2)

def neighbor_align_binned(thetas_dict, neighbor_chart, T,bin_num):
    keys = thetas_dict.keys()
    total_alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
    for key in keys:
        alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
        for neighbor in neighbor_chart[key]:
            alignment += theta_dot(thetas_dict[key], thetas_dict[neighbor])
        total_alignment += alignment
    T_clean= T[1:]
    align_clean = total_alignment
    align_binned= []
    T_binned = []
    k = 0
    align_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            align_accum.append(align_clean[i])
            T_accum.append(T_clean[i])
            k+=1
            continue
        if k == bin_num:
            align_binned.append(np.average(align_accum))
            T_binned.append(np.average(T_accum))
            align_accum = []
            T_accum = []
            k = 0
            continue
    return T_binned, align_binned
    


def hxsq_cnct_align(thetas_dict, neighbor_chart):
    keys = thetas_dict.keys()
    total_alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
    for key in keys:
        alignment = np.zeros(np.shape(thetas_dict[list(keys)[0]])[0])
        for neighbor in neighbor_chart[key]:
            if (neighbor[1] != key[1]) or (neighbor[4] != key[4]):
                alignment += theta_dot(thetas_dict[key], thetas_dict[neighbor])
        total_alignment += alignment
    return total_alignment / (len(list(keys)) * 2) #* 2)

def specific_heat(E, T, half_bin_num):
    bin_num = half_bin_num*2
    Cv = []
    T_clean= []
    E_clean= []
    for i in range(len(E)-1):
        if ((E[i] == None) or (E[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
            continue
        else:
            Cv.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
            T_clean.append(T[i])
            E_clean.append(E[i])
    Cv_floating_avg = []
    T_floating_avg = []
    for i in range(len(Cv) - bin_num):
        Cv_floating_avg.append(np.average(Cv[i-half_bin_num:i +half_bin_num]))
        T_floating_avg.append(np.average(T_clean[i-half_bin_num:i + half_bin_num]))
    E_binned= []
    T_binned = []
    k = 0
    E_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            E_accum.append(E_clean[i])
            T_accum.append(T_clean[i])
            k+=1
        if k ==bin_num:
            E_binned.append(np.average(E_accum))
            T_binned.append(np.average(T_accum))
            E_accum = []
            T_accum = []
            k = 0
    Cv_binned=[]
    for i in range(len(E_binned) - 1):
        Cv_binned.append((E_binned[i] - E_binned[i+1]) / (T_binned[i] - T_binned[i+1]))
    return T_floating_avg, Cv_floating_avg

def specific_heat_binned(E, T, half_bin_num, kb):
    bin_num = half_bin_num*2
    Cv = []
    T_clean= []
    E_clean= []
    for i in range(len(E)-1):
        if ((E[i] == None) or (E[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
            continue
        else:
            Cv.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
            T_clean.append(T[i])
            E_clean.append(E[i])
    E_std_binned= []
    T_binned = []
    k = 0
    E_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            E_accum.append(E_clean[i])
            T_accum.append(T_clean[i])
            k+=1
            continue
        if k == bin_num:
            E_std_binned.append(np.std(E_accum))
            T_binned.append(np.average(T_accum))
            E_accum = []
            T_accum = []
            k = 0
            continue
    Cv_binned=[]
    for i in range(len(E_std_binned) - 1):
        Cv_binned.append((E_std_binned[i]/ T_binned[i])**2 / kb)
    return T_binned[:-1], Cv_binned

def Total_Magnetization(data, T, bin_num):
    all_magnetizations = []
    first = True
    for key in data.keys():
        if first:
            total_vec_x, total_vec_y =  np.cos(data[key]), np.sin(data[key])
            first=False
        else:
            total_vec_x += np.cos(data[key])
            total_vec_y += np.sin(data[key])
    mag_sqr = (total_vec_x**2 + total_vec_y**2)
    mag_clean= mag_sqr
    T_clean= T[1:]

    mag_binned= []
    T_binned = []
    k = 0
    mag_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            mag_accum.append(mag_clean[i])
            T_accum.append(T_clean[i])
            k+=1
            continue
        if k == bin_num:
            mag_binned.append(np.average(mag_accum))
            T_binned.append(np.average(T_accum))
            mag_accum = []
            T_accum = []
            k = 0
            continue
    return T_binned, mag_binned

def Binned_Energy(Energy, T, bin_num):
    T_clean= T[1:]
    E_clean = Energy
    E_binned= []
    T_binned = []
    k = 0
    E_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            E_accum.append(E_clean[i])
            T_accum.append(T_clean[i])
            k+=1
            continue
        if k == bin_num:
            E_binned.append(np.average(E_accum))
            T_binned.append(np.average(T_accum))
            E_accum = []
            T_accum = []
            k = 0
            continue
    return T_binned, E_binned
def susceptibility(data, T, bin_num):
    all_magnetizations = []
    first = True
    for key in data.keys():
        if first:
            total_vec_x, total_vec_y =  np.cos(data[key]), np.sin(data[key])
            first=False
        else:
            total_vec_x += np.cos(data[key])
            total_vec_y += np.sin(data[key])
    mag_clean = total_vec_x
    T_clean= T[1:]

    chi_binned= []
    T_binned = []
    k = 0
    mag_accum = []
    T_accum = []
    for i in range(len(T_clean)):
        if k < bin_num:
            mag_accum.append(mag_clean[i])
            T_accum.append(T_clean[i])
            k+=1
            continue
        if k == bin_num:
            chi_binned.append((np.std(mag_accum)**2) / np.average(T_accum))
            T_binned.append(np.average(T_accum))
            mag_accum = []
            T_accum = []
            k = 0
            continue
    #print(T_binned, mag_binned)
    return T_binned, chi_binned
# def Mag_binned(data, T, half_bin_num):
#     bin_num = half_bin_num*2
#     Mag = []
#     T_clean= []
#     for i in range(len(data)-1):
#         if ((data[i] == None) or (data[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
#             continue
#         else:
#             Mag.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
#             T_clean.append(T[i])
#     E_std_binned= []
#     T_binned = []
#     k = 0
#     E_accum = []
#     T_accum = []
#     for i in range(len(T_clean)):
#         if k < bin_num:
#             E_accum.append(E_clean[i])
#             T_accum.append(T_clean[i])
#             k+=1
#             continue
#         if k == bin_num:
#             E_std_binned.append(np.std(E_accum))
#             T_binned.append(np.average(T_accum))
#             E_accum = []
#             T_accum = []
#             k = 0
#             continue
#     Cv_binned=[]
#     for i in range(len(E_std_binned) - 1):
#         Cv_binned.append((E_std_binned[i]/ T_binned[i])**2)
#     return T_binned[:-1], Cv_binned

def specific_heat_binned_by_T(E, T, T_step):
    Cv = []
    T_clean= []
    E_clean= []
    for i in range(len(E)-1):
        if ((E[i] == None) or (E[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
            continue
        else:
            Cv.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
            T_clean.append(T[i])
            E_clean.append(E[i])
    E_std_binned= []
    T_binned = []
    E_accum = []
    T_accum = []
    T_step_so_far = 0
    for i in range(len(T_clean) - 1):
        if T_step_so_far < T_step:
            E_accum.append(E_clean[i])
            T_accum.append(T_clean[i])
            T_step_so_far += (T_clean[i] - T_clean[i+1])
            continue
        if T_step_so_far >= T_step:
            E_std_binned.append(np.std(E_accum))
            T_binned.append(np.average(T_accum))
            E_accum = []
            T_accum = []
            T_step_so_far = 0
            continue
    Cv_binned=[]
    for i in range(len(E_std_binned) - 1):
        Cv_binned.append((E_std_binned[i] / T_binned[i])**2)
    return T_binned[:-1], Cv_binned



def specific_heat_binned_by_T_old(E, T, T_step):
    Cv = []
    T_clean= []
    E_clean= []
    for i in range(len(E)-1):
        if ((E[i] == None) or (E[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
            continue
        else:
            Cv.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
            T_clean.append(T[i])
            E_clean.append(E[i])
    E_binned= []
    T_binned = []
    E_accum = []
    T_accum = []
    T_step_so_far = 0
    for i in range(len(T_clean) - 1):
        if T_step_so_far < T_step:
            E_accum.append(E_clean[i])
            T_accum.append(T_clean[i])
            T_step_so_far += (T_clean[i] - T_clean[i+1])
            continue
        if T_step_so_far >= T_step:
            E_binned.append(np.average(E_accum))
            T_binned.append(np.average(T_accum))
            E_accum = []
            T_accum = []
            T_step_so_far = 0
            continue
    Cv_binned=[]
    for i in range(len(E_binned) - 1):
        Cv_binned.append((E_binned[i] - E_binned[i+1]) / (T_binned[i] - T_binned[i+1]))
    return T_binned[:-1], Cv_binned
def specific_heat_2(E, T, half_bin_num):
    bin_num = half_bin_num*2
    Cv = []
    T_clean= []
    E_clean= []
    for i in range(len(E)-1):
        if ((E[i] == None) or (E[i+1] == None)) or ((T[i] == None) or (T[i+1] == None)):
            continue
        else:
            Cv.append((E[i] - E[i+1]) / (T[i]- T[i+1]))
            T_clean.append(T[i])
            E_clean.append(E[i])
    T_floating_avg = []
    E_floating_avg = []
    for i in range(len(T) - bin_num):
        E_floating_avg.append(np.average(E_clean[i-half_bin_num:i + half_bin_num]))
        T_floating_avg.append(np.average(T_clean[i-half_bin_num:i + half_bin_num]))
    Cv_for_E_floating_avg = []
    for i in range(len(T_floating_avg) -1):
        Cv_for_E_floating_avg.append((E_floating_avg[i] - E_floating_avg[i+1]) / (T_floating_avg[i]- T_floating_avg[i+1]))

    return T_floating_avg[:-1], Cv_for_E_floating_avg

# def upscore3(thetas):
#     vect = lambda thet: np.array([np.sin(thet), np.cos(thet)]) 
#     vectors = np.transpose(np.array([vect(theta) for theta in thetas]))
#     score = np.zeros(np.shape(vectors)[0])
#     def dot2(m4):
#         m22 = np.reshape(m4, (2,2))
#         return np.dot(m22[0], m22[1])
#     for i in range(6):
#         conc = np.concatenate([vectors[:,:,i], vectors[:,:,(i + 1)%6]], axis=1)
#         new_score = np.apply_along_axis(dot2, 1, conc)
#         score += new_score
#     return score / 6
