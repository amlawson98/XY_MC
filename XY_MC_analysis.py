from tkinter import *
import numpy as np 
import random
import copy
from matplotlib import pyplot as plt
from math import pi as pi
import math
from PIL import ImageGrab
import os
import pickle
import tqdm
from alignment_funcs import susceptibility, neighbor_align_binned, Binned_Energy, hex_align, Total_Magnetization, neighbor_align, hxsq_cnct_align, specific_heat, specific_heat_2, specific_heat_binned, specific_heat_binned_by_T
from particle_utils import string_to_array, arr_of_arrs_to_arr_of_arrs_of_arrs

plt.rcParams.update({'font.size': 16})
x_max = 60
save_nametr = 'triangular_lattice_64-pts_100000-epochs'
save_namehxsqtr = 'hxsqtr_lattice_96-pts_100000-epochs'
#save_namekag = 'kagome_lattice_108-pts_100000-epochs'
save_namekag = 'kagome_lattice_48-pts_100000-epochs_120-ticks_2'
save_name_kag_2 = 'kagome_lattice_108-pts_1-J_1000-epochs_120-ticks_25-Tmax_0.001-Tmin'
save_name_snub_trhx = 'snubtrhx_lattice_96-pts_100000-epochs'
save_name_snub_sqr = 'snub_square_lattice_128-pts_100000-epochs_120-ticks'
save_name_hex = 'hexagonal_lattice_36-pts_100000-epochs_120-ticks'
save_name_square = 'square_lattice_64-pts_100000-epochs_120-ticks'
#save_name_square_high = 'square_lattice_64-pts_1-J_1000-epochs_120-ticks_25-Tmax_1-Tmin'

#save_name_square ='square_lattice_16-pts_1000-epochs_120-ticks_100-Tmax'
#save_name_square = 'square_lattice_16-pts_-1-J_10000-epochs_120-ticks_200-Tmax'
save_name_trunc = 'trunc_hex_lattice_216-pts_100000-epochs_120-ticks_2'
# with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
#     data, neighbor_chart, T, E = pickle.load(f)
#     N = len(neighbor_chart.keys())

def comb_data(max_temp, data, T, E):
    for i in range(len(T)):
        if T[i] is not None:
            if T[i] > max_temp:
                for k in data.keys():
                    data[k][i] = None
                E[i] = None
                T[i] = None
        else:
            for k in data.keys():
                data[k][i] = None
            E[i] = None
            T[i] = None
    for k in data.keys():
        data[k] = [dki for dki in data[k] if dki is not None]
    T = [t for t in T if t is not None]
    E = [e for e in E if e is not None]
    return data, T, E




def composite_angle(six_angles):
    base_vec = np.array([0.,0.])
    k = 0
    for angle in six_angles:
        if k in [0, 2, 4]:
            base_vec += np.array([np.cos(angle), np.sin(angle)])
        else:
            base_vec -= np.array([np.cos(angle), np.sin(angle)])
        k+=1
    main_vec = base_vec / 6
    thet = np.arctan2(base_vec[1], base_vec[0])
    return thet


def hxsqtr_to_triang_lattice(hxsqrtr_data, hxsqtr_neighbors):
    arr_keys = []
    for key in hxsqtr_neighbors.keys():
        arr_keys.append(string_to_array(key))
    arr_of_hexagon_indices, dict_of_hexagon_indices = arr_of_arrs_to_arr_of_arrs_of_arrs(arr_keys)
    new_neighbor_dict = {}
    for key in dict_of_hexagon_indices.keys():
        if key in new_neighbor_dict.keys():
            for site in dict_of_hexagon_indices[key]:
                for neighbor in hxsqtr_neighbors[str(site)]:
                    if str(string_to_array(neighbor)[:-1]) not in new_neighbor_dict[key] and ((str(string_to_array(neighbor)[:-1])) != key):
                        new_neighbor_dict[key].append(str(string_to_array(neighbor)[:-1]))
        else:
            all_base_neighbors = []
            for site in dict_of_hexagon_indices[key]:
                for neighbor in hxsqtr_neighbors[str(site)]:
                    if [(str(string_to_array(neighbor)[:-1]))] not in all_base_neighbors and ((str(string_to_array(neighbor)[:-1])) != key):
                        all_base_neighbors.append((str(string_to_array(neighbor)[:-1])))
            new_neighbor_dict[key] = all_base_neighbors
    
    new_data = {}
    for key in dict_of_hexagon_indices.keys():
        base_spins = np.zeros((6, np.shape(hxsqrtr_data['[0, 0, 0]'])[0]))
        for k in range(len(dict_of_hexagon_indices[key])):
            base_spins[k] = hxsqrtr_data[str(dict_of_hexagon_indices[key][k])]
        composite_spins = np.apply_along_axis(composite_angle, 1, np.transpose(base_spins))
        new_data[key] = composite_spins
    return new_data, new_neighbor_dict



def plot_energy():
    plt.title(save_name)
    plt.plot(np.arange(len(E)), E)
    plt.xlabel('number of iterations')
    plt.ylabel('Energy')
    #plt.show()

def divide_nones(x1, x2):
    if x1 == None or x2 == None:
        return None
    else:
        return x1/x2

def plot_S(save_name, num, kb):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(200, data, T, E)
    T_arr, Cv_arr = specific_heat_binned(E, np.array(T), num, kb)
    #Cv = np.vectorize(divide_nones)(np.array(E), np.array(T))

    T_mid_points = [(T_arr[i] + T_arr[i + 1])/2 for i in range(len(T_arr[:-1]))]
    T_diffs = [(T_arr[i] - T_arr[i + 1]) for i in range(len(T_arr[:-1]))]
    Cv_mid_points= [(Cv_arr[i] + Cv_arr[i + 1])/2 for i in range(len(Cv_arr[:-1]))]

    dS_s = np.array(Cv_mid_points) * T_diffs/ np.array(T_mid_points)

    S = -np.array([np.sum(dS_s[:i]) for i in range(len(dS_s))])

    plt.plot(np.array(T_mid_points)*kb, S / N, label =save_name)
    #plt.scatter(T_binned[:-1], Cv_binned, label ="Cv", alpha = .1)
    plt.xlim(0, x_max)
    plt.xlabel(r'Temperature ($k_b = 1$)')
    plt.ylabel(r'$\frac{\Delta S}{N}$ Entropy per Site')

def plot_Cv_T(save_name, num, kb):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(200, data, T, E)
    T_arr, Cv_arr = specific_heat_binned(E, np.array(T), num, kb)
    #Cv = np.vectorize(divide_nones)(np.array(E), np.array(T))
    plt.plot(np.array(T_arr)*kb, np.array(Cv_arr) / N, label =save_name)
    #plt.scatter(T_binned[:-1], Cv_binned, label ="Cv", alpha = .1)
    plt.xlim(0, x_max)
    plt.xlabel(r'Temperature ($k_b = 1$)')
    plt.ylabel(r'$C_v$ / $N$ Specific Heat per Site')

def plot_magnetization(save_name, bin_num):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(15, data, T, E)
    T_arr, M_arr = Total_Magnetization(data, np.array(T), bin_num)
    plt.plot(T_arr, np.array(M_arr) / N, label =save_name)
    plt.xlim(0, x_max)
    plt.xlabel('temperature*kb')
    plt.ylabel('<M^2>/N')

def plot_susceptibility(save_name, bin_num):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(15, data, T, E)
    T_arr, chi_arr = susceptibility(data, np.array(T), bin_num)
    plt.plot(T_arr, (np.array(chi_arr)) / N, label =save_name)
    plt.xlim(0, x_max)
    plt.xlabel('temperature*kb')
    plt.ylabel('(Susceptibility)/N')

def plot_hex_align_E():
    plt.title(save_name)
    plt.scatter(np.array(E)/N, hex_align(data))
    plt.xlabel('Energy/N')
    plt.ylabel('hexagonal order')
    #plt.show()

def plot_hex_align_T():
    plt.title(save_name)
    plt.scatter(T, hex_align(data), label="hxsq_align")
    plt.xlabel('temperature')
    plt.ylabel('hexagonal order')
    #plt.show()

def plot_neighbor_align_T(save_name, half_bin_num):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(15, data, T, E)
    T, neighbor_align = neighbor_align_binned(data, neighbor_chart, T, half_bin_num)
    plt.xlabel('temperature')
    plt.ylabel('neighbor alignment/ N')
    plt.plot(T, np.array(neighbor_align) / N, label =save_name)

def plot_hxsq_cnct_align_T():
    plt.title(save_name)
    plt.scatter(T, hxsq_cnct_align(data, neighbor_chart), label="hxsq_cnct_align")
    plt.xlabel('temperature*kb')
    plt.ylabel('neighbor order')
    #plt.show()

def plot_Energy(save_name, bin_num):
    with open('lattice_data/' + save_name + '.pkl', 'rb') as f:
        data, neighbor_chart, T, E = pickle.load(f)
    N = len(neighbor_chart.keys())
    data, T, E = comb_data(15, data, T, E)
    T_arr, E_arr = Binned_Energy(E, np.array(T), bin_num)
    print(E_arr[-1]/N)
    plt.plot(T_arr, np.array(E_arr) / N, label =save_name)
    plt.xlim(0, x_max)
    plt.xlabel(r'Temperature ($k_b = 1$)')
    plt.ylabel(r'$E$ / $N$ Energy Heat per Site')

#plot_neighbor_align_T(data, neighbor_chart, T)
#plt.show()
#plot_hex_align_T()
#plot_hxsq_cnct_align_T()
plt.xscale('log')
bin_num = 1000
kb = 1
#plot_Cv_T(save_nametr, bin_num, alpha=1)
#plot_Cv_T(save_namehxsqtr, bin_num, alpha=1)

# plot_Cv_T(save_name_snub_trhx, bin_num, kb)
# plot_Cv_T(save_name_snub_sqr, bin_num, kb)
# plot_Cv_T(save_name_hex, bin_num, kb)

plot_S(save_namekag, bin_num, kb)
plot_S(save_name_square, bin_num, kb)
plot_S(save_name_trunc, bin_num, kb)


# plot_Energy(save_namekag, bin_num)
# plot_Energy(save_name_square, bin_num)
# plot_Energy(save_name_trunc, bin_num)

J=1
T_c = abs(J)*pi/(2*kb)
plt.title('Entropy vs. Temperature for Kagome, Truncated Hexagonal, and Square Anti-Ferromagnets')
plt.legend()
plt.show()

# plot_Cv_T(save_name_trunc, bin_num, alpha=1)
# # comp_tri_data, comp_tri_neighbors = hxsqtr_to_triang_lattice(data, neighbor_chart)

# # plot_neighbor_align_T(comp_tri_data, comp_tri_neighbors, T)

# plt.title("Specific Heat vs. Temperature for 4 lattices")
# plt.xscale('log')
# plot_susceptibility(save_name_square, bin_num)
# plot_susceptibility(save_namekag, bin_num)
# plot_susceptibility(save_name_trunc, bin_num)
# plt.legend()
# plt.show()

# bin_num = 300
# plt.xscale('log')
# plot_Energy(save_name_square, bin_num)
# plot_Energy(save_namekag, bin_num)
# plot_Energy(save_name_trunc, bin_num)
# plt.legend()
# plt.show()

# bin_num = 300
# plt.xscale('log')
# plot_neighbor_align_T(save_name_square, bin_num)
# plot_neighbor_align_T(save_namekag, bin_num)
# plot_neighbor_align_T(save_name_trunc, bin_num)
# plt.legend()
# plt.show()
# # color = plt.cm.viridis()
# ax4.scatter(x =np.array(y),y =np.array(hx)  / np.array(y), c=x, cmap = 'viridis')
# ax4.plot(np.array(y), np.array(hx)  / np.array(y), linewidth = .5, alpha = .5)
# ax4.set_ylabel('hexagonal disorder', color='r')
# ax4.set_xlabel('Energy', color='b')
# #ax4.set_zlabel('time (epochs)')
# # plt.show()
# print(np.shape(hx), np.shape(ng))
# ax4.scatter(x =np.array(y),y =np.array(hx), c=x, cmap = 'viridis')
# ax4.plot(np.array(y), np.array(hx), linewidth = .5, alpha = .5)
# ax4.set_ylabel('hexagonal disorder', color='r')
# ax4.set_xlabel('Energy', color='b')
# #ax4.set_zlabel('time (epochs)')
# plt.show()

# fig3 = plt.figure()
# ax5 = fig3.add_subplot(111)
# ax5.scatter(x =np.array(y), y =np.array(ng), c=x, cmap = 'viridis')
# ax5.plot(np.array(y), np.array(ng), linewidth = .5, alpha = .5)
# ax5.set_ylabel('neighbor_disorder', color='r')
# ax5.set_xlabel('Energy', color='b')
# plt.show()

# fig4 = plt.figure()
# ax6 = fig4.add_subplot(111)
# ax6.scatter(x =np.array(y),y =np.array(hx) - np.array(ng), c=x, cmap = 'viridis')
# ax6.plot(np.array(y), np.array(hx) - np.array(ng), linewidth = .5, alpha = .5)
# ax6.set_ylabel('hex_disorder - neighbor_disorder', color='r')
# ax6.set_xlabel('Energy', color='b')
# plt.show()
# closed = False
# if not interact:
#     plt.close('all')
# if save_best:
#     try:
#         with open('lattice_dicts/' + save_name + '.pkl', 'rb') as f:
#             old_best = pickle.load(f)
#             if Total_Energy(old_best, Energy) >= Total_Energy(best_config, Energy):
#                 better = True
#                 print("new best Energy: ", Total_Energy(best_config, Energy))             
#                 with open('lattice_dicts/'+ save_name + '.pkl', 'wb') as f:
#                     pickle.dump(best_config, f, pickle.HIGHEST_PROTOCOL)
#             else:
#                 better = False
#     except:
#         better = True
#         print("first best Energy: ", Total_Energy(best_config, Energy))             
#         with open('lattice_dicts/'+ save_name + '.pkl', 'wb') as f:
#             pickle.dump(best_config, f, pickle.HIGHEST_PROTOCOL)

