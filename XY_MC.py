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
from numba import jit
from particle_utils import draw_arrow, size, string_to_array
from lattices import hxsqtr_lattice, triangular_lattice, kagome_lattice, snbtrhx_lattice, hexagonal_lattice, trunc_hex_lattice, square_lattice

#cmap = plt.get_cmap('gist_rainbow')
window_size = 700.0 

class spin(object):
    def __init__(self, name,  neighbors, pos, direction= random.uniform(0,2*pi), magnitude = 60):
        self.name = name
        self.direction = direction
        self.neighbors = neighbors
        self.pos = pos

        self.color = 'light blue'
        self.magnitude = magnitude

    def show(self):
        print(self.name, self.direction, self.neighbors, self.pos)

    def radius(self):
        return 10

    def get_corners(self):
        xy0 = np.array(
            [self.pos[0] - self.radius(), self.pos[1] + self.radius()])
        xy1 = np.array(
            [self.pos[0] + self.radius(), self.pos[1] - self.radius()])
        return xy0, xy1

    def mag_moment_xy(self):
        return self.magnitude * np.array([np.cos(self.direction), np.sin(self.direction)])

    def draw(self, canvas):
        xy0, xy1 = self.get_corners()
        # k = self.name[7]
        # if k == "0":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
        #                    outline="black", fill='orange')
        #     canvas.pack()
        # if k == "1":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
        #                    outline="black", fill='green')
        #     canvas.pack()
        # if k == "2":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
        #                    outline="black", fill='red')
        #     canvas.pack()
        # if k == "3":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
        #                    outline="black", fill='light blue')
        #     canvas.pack()
        # if k == "4":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
        #                    outline="black", fill='purple')
        #     canvas.pack()
        # if k == "5":
        #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1], outline = "black", fill='pink')
        #     canvas.pack()
        draw_arrow(canvas, self.pos - (self.mag_moment_xy() / 2), self.mag_moment_xy())

    def draw_lines(self, spins_dict):
        for name in self.neighbors:
            if size(self.pos - spins_dict[name].pos) < (window_size / 2):
                midpoint = self.pos + (spins_dict[name].pos - self.pos) / 2
                canvas.create_line(self.pos[0], self.pos[1], midpoint[0], midpoint[1], fill= 'blue')


    def energy(self, others, Hamiltonian_particle):
        return Hamiltonian_particle(self, others)

def display_grid(number, grid_color):
    canvas.delete("all")
    for x in np.arange(window_size / (2 * number), window_size + (window_size /(2 *  number)), window_size / number):
        canvas.create_line(x, 0, x, window_size, fill=grid_color)
    for y in np.arange(window_size / (2 * number), window_size + (window_size / (2 * number)), window_size / number):
        canvas.create_line(0, y, window_size, y, fill=grid_color)

# def grid(number = 4, sig_omega = 0):
#     particles = np.empty((number, number), spin)
#     for i in range(number):
#         for j in range(number):
#             particles[i, j] = spin(
#                 pos = np.array([(window_size / (2 * number)) + i * window_size / number, (window_size / (2 * number)) + j * window_size / number]),
#                 name = [i, j],
#                 direction = random.uniform(0, 2 * pi),
#                 omega = random.uniform(-1 * sig_omega, sig_omega))
#     global glob_num_in_line
#     glob_num_in_line = number 
#     return particles.reshape(number ** 2)


def display_during(spins_dict, num_lines, canvas, save_name = "dummy_name", time = 10000, draw_lines = True, interact = True, former_best_energy=0):
    canvas.pack()
    if Total_Energy(spins_dict, Energy) < former_best_energy:
        for j in range(5):
        #display_grid(num_lines, "gray")
            for k in spins_dict.keys():
                spins_dict[k].draw(canvas)
                if draw_lines:
                    spins_dict[k].draw_lines(spins_dict)
            canvas.pack()
            canvas.update()
            if (j >= 0 and j < 5):
                #root.overrideredirect(1)
                #root.geometry('%dx%d+%d+%d' % (window_size, window_size, 0, 0))
                #canvas.pack(side=RIGHT, expand=YES, fill=BOTH)
                canvas.update()
            elif (j == 5):
                print("saving to:", os.getcwd() + "/XY_sim_images")
                savename = str("XY_sim_images/" + save_name).format(1)
                ImageGrab.grab((0,0,window_size,window_size)).save(savename + '.jpg')
                #root.overrideredirect(0)
    else:
        for k in spins_dict.keys():
            spins_dict[k].draw(canvas)
            if draw_lines:
                spins_dict[k].draw_lines(spins_dict)
            canvas.pack()
            canvas.update()

def display_last(spins_dict, num_lines, canvas, save_best = True, save_name = "dummy_name", time = 10000, draw_lines = True, interact = True):
    canvas.pack()
    for j in range(time):
        #display_grid(num_lines, "gray")
        for k in spins_dict.keys():
            spins_dict[k].draw(canvas)
            if draw_lines:
                spins_dict[k].draw_lines(spins_dict)
        canvas.pack()
        canvas.update()
        if save_best and (j >= 0 and j < 5):
            #root.overrideredirect(1)
            #root.geometry('%dx%d+%d+%d' % (window_size, window_size, 0, 0))
            #canvas.pack(side=RIGHT, expand=NO, fill=NO)
            canvas.update()
        elif save_best and (j == 5):
            print("saving to:", os.getcwd() + "/XY_sim_images")
            savename = str("XY_sim_images/" + save_name).format(1)
            #ImageGrab.grab((0,0,window_size,window_size)).save(savename + '.jpg')
            #root.overrideredirect(0)
        elif not interact:
            root.destroy()
            break

def Energy(main_spin, spins_dict, J = 1, H = 0, main_theta_shift = 0):
    directions = list(map(lambda name: spins_dict[name].direction, main_spin.neighbors))
    main_theta = main_spin.direction + main_theta_shift
    #main_spin.show()
    #print(main_theta)
    #print(directions)
    E = H * np.cos(main_theta)
    for theta in directions:
        E += (J * np.cos(main_theta - theta))
        #print(E)
    return E

def Total_Energy(spins_dict, Energy_func):
    E_tot = 0
    for spin_key in spins_dict.keys():
        E_tot += Energy_func(spins_dict[spin_key], spins_dict)
    return E_tot

def Best_Random_Config(spins_dict, generate_spins, Energy_func, num_tries_left):
    if num_tries_left == 0:
        return spins_dict
    else:
        new_spins = generate_spins()
        #print(Total_Energy(spins_dict, Energy_func))
        #print(Total_Energy(new_spins, Energy_func))
        if Total_Energy(spins_dict[0], Energy_func) <= Total_Energy(new_spins[0], Energy_func):
            best_spins = copy.copy(spins_dict)
        else:
            best_spins = copy.copy(new_spins)
        #print(Total_Energy(best_spins, Energy_func))
        return Best_Random_Config(best_spins, generate_spins, Energy_func, num_tries_left - 1)

num_lines = 8

chance_min = .8
amount_flipped = .5

@jit
def upsate_spin_step(spin, spin_dict, Energy_func, beta, num_ticks):
    chance = min(chance_min + ((1 - chance_min) * 1.5 * (i / epochs)), 1)
    cur_Energy = Energy_func(cur_spin, spins_dict)
    #(np.log(i + 100))
    random_shift = (random.uniform(0, 1)**.5)
    plus_spin = copy.copy(cur_spin)
    plus_spin.direction += update_step * random_shift
    plus_Energy = Energy_func(plus_spin, spins_dict)

    minus_spin = copy.copy(cur_spin)
    minus_spin.direction -= update_step * random_shift
    minus_Energy = Energy_func(minus_spin, spins_dict)

    plus_best = (plus_Energy < cur_Energy) and (plus_Energy < minus_Energy)   
    cur_best = (cur_Energy < plus_Energy) and (cur_Energy < minus_Energy)  
    minus_best = (minus_Energy < plus_Energy) and (cur_Energy > minus_Energy) 

    if plus_best:
        if random.uniform(0,1) < chance:
            cur_spin = plus_spin
        else:
            if (minus_Energy < cur_Energy):
                if random.uniform(0,1) < chance:
                    cur_spin = minus_spin
                else:
                    cur_spin = cur_spin
            else:
                if random.uniform(0,1) < chance:
                    cur_spin = cur_spin
                else:
                    cur_spin = minus_spin

    if minus_best:
        if random.uniform(0,1) < chance:
            cur_spin = minus_spin
        else:
            if (plus_Energy < cur_Energy):
                if random.uniform(0,1) < chance:
                    cur_spin = plus_spin
                else:
                    cur_spin = cur_spin
            else:
                if random.uniform(0,1) < chance:
                    cur_spin = cur_spin
                else:
                    cur_spin = plus_spin

    if cur_best:
        if random.uniform(0,1) < chance:
            cur_spin = cur_spin
        else:
            if (plus_Energy < minus_Energy):
                if random.uniform(0,1) < chance:
                    cur_spin = plus_spin
                else:
                    cur_spin = minus_spin
            else:
                if random.uniform(0,1) < chance:
                    cur_spin = minus_spin
                else:
                    cur_spin = plus_spin
    spins_dict[spin_key] = cur_spin

linger = 0
linger_temp = 0.05
temp_exponent = 1
def temperature(epochs, epochs_total):
    if ((epochs_total - epochs) > linger):
        epoch_effective = epochs * (epochs_total/ (epochs_total - linger)) + 0.00000001
        #print(epochs, epoch_effective)
        return (((epochs_total/(epoch_effective)) - 1)) ** temp_exponent
    else:
        return linger_temp

tau = .1
T_0 = 10
def exp_temp(epochs, epochs_total):
    return T_0 * np.exp(-1*(epochs/ epochs_total) / tau)

kb = .1
def update_spin_statistical(spin, spin_dict, Energy_func, temp, num_ticks):
    possible_thetas = []
    possible_energies = []
    for j in range(num_ticks):
        main_theta = spin.direction
        j_theta_shift = (j / num_ticks)*2*pi
        possible_thetas.append(j_theta_shift + main_theta)
        possible_energies.append(Energy_func(spin, spin_dict, j_theta_shift))
    relative_probabilities = np.exp(-np.array(possible_energies)/(temp * kb))
    Z = np.sum(relative_probabilities)
    #print(relative_probabilities)
    #print(temp)
    probabilities = relative_probabilities / Z
    #print(probabilities)
    prob_num = random.uniform(0,1)
    #print(prob_num)
    index = 0
    for j in range(num_ticks):
        if prob_num < probabilities[j]:
            index = j
            break
        else:
            prob_num -= probabilities[j]
    new_theta = possible_thetas[index]
    spin.direction = new_theta

@jit
def sweep_update(spins_dict, Energy_func, num_sweeps_left, cur_epochs, total_epochs, 
    update_step = 1, i = 1, show = False, former_best_energy = 0, num_ticks =4):
    if num_sweeps_left == 0:
        if show:
            canvas.delete("all")
            display_during(spins_dict, num_lines, canvas, former_best_energy = former_best_energy) 
            canvas.pack()
            root.update()
        return spins_dict
    else:
        temp =exp_temp(cur_epochs, total_epochs)
        all_keys = spins_dict.keys()
        some_keys = all_keys
        #random.sample(all_keys, math.floor(((num_lines ** 2) * amount_flipped)))
        for spin_key in some_keys:
            cur_spin = spins_dict[spin_key]
            update_spin_statistical(cur_spin, spins_dict, Energy_func, temp, num_ticks)
        return sweep_update(spins_dict, Energy_func, num_sweeps_left - 1, cur_epochs, total_epochs,
         update_step, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks)

J = -1
H = 0

# def upscore(top_thet, top_left_thet, bottom_left_thet, bottom_thet, bottom_right_thet, top_right_thet):
#     return np.var([top_thet%(2*pi), 
#         (pi + top_left_thet) % (2*pi),
#         bottom_left_thet % (2*pi),
#         (pi + bottom_thet) %(2*pi),
#         bottom_right_thet %(2*pi),
#         (pi + top_right_thet) %(2*pi)
#      ]) / (5*(pi**2)/4)

# def upscore2(top_thet, top_left_thet, bottom_left_thet, bottom_thet, bottom_right_thet, top_right_thet):
#     thetas = np.array([top_thet, top_left_thet, bottom_left_thet, bottom_thet, bottom_right_thet, top_right_thet])
#     sins_coss = np.array([np.sin(thetas), np.cos(thetas)])
#     vectors = [np.array([np.sin(theta), np.cos(theta)]) for theta in thetas]
#     score = 0
#     for i in range(6):
#         score += np.dot(vectors[i], vectors[(i + 1)%6])
#     return score / 6
# pi =np.pi

# def alignment_of_hexagon(hex_spins):
#     top_thet = hex_spins[0].direction
#     top_left_thet = hex_spins[1].direction
#     bottom_left_thet = hex_spins[2].direction
#     bottom_thet = hex_spins[3].direction
#     bottom_right_thet = hex_spins[4].direction
#     top_right_thet = hex_spins[5].direction

#     #print(+++++++++++)
#     #print(top_thet, top_left_thet, top_right_thet, bottom_left_thet, bottom_right_thet)
#     return upscore2(top_thet, top_left_thet, bottom_left_thet, bottom_thet, bottom_right_thet, top_right_thet)

# def hex_align_of_hxsqrtr_lattice(hxsqrtr_dict):
#     #print(hxsqrtr_dict.keys())
#     keys = hxsqrtr_dict.keys()
#     keys_as_arrays = [string_to_array(key) for key in keys] 
#     keys_no_ending = np.unique([arr_key[:-1] for arr_key in keys_as_arrays], axis = 0)
#     #print(keys_no_ending)
#     total_alignment = 0
#     for block in keys_no_ending:
#         #print(str(np.append(block, [4])))
#         block_alignment = alignment_of_hexagon([hxsqrtr_dict[str(list(np.append(block, [ending])))] for ending in range(6)])
#         total_alignment += block_alignment
#     return total_alignment / len(keys_no_ending)

def alignment_of_spin(key, lattice_dict):
    main_spin = lattice_dict[key]
    neighbors = main_spin.neighbors
    main_theta = lattice_dict[key].direction
    alignment = 0
    for neighbor in neighbors:
        neighbor_theta = lattice_dict[neighbor].direction
        alignment += np.dot(np.array([np.cos(main_theta), np.sin(main_theta)]), np.array([np.cos(neighbor_theta), np.sin(neighbor_theta)]))
    return alignment / len(neighbors)


def neighbor_align(lattice_dict):
    keys = lattice_dict.keys()
    total_alignment =0
    for key in keys:
        total_alignment += alignment_of_spin(key, lattice_dict)
    return total_alignment / (len(keys))

def init_data(lattice_dict):
    keys = lattice_dict.keys()
    data = {}
    neighbor_chart = {}
    for key in keys:
        data[key] = np.array([copy.copy(lattice_dict[key].direction)])
        neighbor_chart[key] = copy.copy(lattice_dict[key].neighbors)
    return data, neighbor_chart

def update_data(old_data, new_lattice):
    for key in new_lattice.keys():
        old_data[key] = np.append(old_data[key], copy.copy(new_lattice[key].direction))
    return old_data

def save_data(data, neighbor_chart, T, E, place, name):
    with open(place + name + '.pkl', 'wb') as f: 
        pickle.dump((data, neighbor_chart, T, E), f, pickle.HIGHEST_PROTOCOL)

sweeps_per_epoch = 10


def optimize_and_draw(lattice, num_lines, num_ticks = 4, save_best = True, save_name = "dummy_name",
 epochs = 5000, interact = True, show = False):
    best_init_config, lattice_type = Best_Random_Config(lattice(num_lines), 
        lambda: lattice(num_lines),
        lambda spin, spins_dict: Energy(spin, spins_dict, J = J, H = H),
        1)
    # with open('lattice_dicts/' + save_name + '.pkl', 'rb') as f:
    #             best_init_config = pickle.load(f)
    # for key in best_init_config.keys():
    #     best_init_config[key].show()
    global root
    root = Tk()
    global canvas
    canvas = Canvas(root, width = window_size, height= window_size)
    best_config = sweep_update(best_init_config, 
        lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift), 
        num_sweeps_left =1, cur_epochs= 1, total_epochs = 2, show=show, num_ticks = num_ticks)
    if save_best:
        try:
            with open('lattice_dicts/' + save_name + '.pkl', 'rb') as f:
                old_best = pickle.load(f)
                former_best_energy = Total_Energy(old_best, Energy)
                if former_best_energy >= Total_Energy(best_config, Energy):
                    better = True
                    print("new best Energy: ", Total_Energy(best_config, Energy))             
                    with open('lattice_dicts/'+ save_name + '.pkl', 'wb') as f:
                        pickle.dump(best_config, f, pickle.HIGHEST_PROTOCOL)
                else:
                    better = False
        except:
            better = True
            former_best_energy = 0
            print("first best Energy: ", Total_Energy(best_config, Energy))             
            with open('lattice_dicts/'+ save_name + '.pkl', 'wb') as f:
                pickle.dump(best_config, f, pickle.HIGHEST_PROTOCOL)
    i = 0
    x = [0]
    y = [Total_Energy(best_init_config, Energy)]
    #ng = [neighbor_align(best_config)]
    # hx = [hex_align_of_hxsqrtr_lattice(best_config)]

    T = [None]
    data, neighbor_chart = init_data(best_init_config)

    global temp
    temp = 1000
    global Eng
    Eng = y[0]
    t = tqdm.trange(epochs, desc = "Start")
    for i in t:
        best_config =  sweep_update(best_init_config, 
        lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift),
        num_sweeps_left = sweeps_per_epoch, cur_epochs = i + 1, total_epochs = epochs, show =show, 
        former_best_energy = former_best_energy, num_ticks = num_ticks)
        Eng = Total_Energy(best_config, lambda spin, spins_dict: Energy(spin, spins_dict, J = J))
        # alignment_hex = hex_align_of_hxsqrtr_lattice(best_config)
        #alignment_ng = neighbor_align(best_config)
        if Eng < former_best_energy:
            better = True
        else:
            better = False
        #display_during(best_config, num_lines, canvas, save = better, save_name = save_name, draw_lines = True, interact = interact)
        y.append(Eng)
        x.append(i)
        # hx.append(alignment_hex)
        #ng.append(alignment_ng)
        temp = exp_temp(i, epochs)
        t.set_description("T = %s, E = %s"%("{:0<14f}".format(temp), "{:0<8f}".format(Eng)))
        data = update_data(data, best_config)
        T.append(temp)  
        #print("T = %s, E = %s"%(round(temp, 3), round(Eng, 3)))

    save_data(data, neighbor_chart, T, y, 'lattice_data/', '%s_%s-pts_%s-epochs'%(lattice_type, len(best_config.keys()), epochs))
    # plt.ion()
    # fig, ax1 = plt.subplots()
    # ax1.plot(x, y, 'b')
    # ax1.set_xlabel('number of iterations')
    # ax1.set_ylabel('Energy', color='b')
    # ax1.tick_params('y', colors='b')

    # ax2 = ax1.twinx()
    # ax2.plot(x, hx, 'r')
    # ax2.set_ylabel('hexagonal disorder', color='r')
    # ax2.tick_params('r', colors='r')

    # ax3 = ax1.twinx()
    # ax3.plot(x, T, 'y')
    # ax3.set_ylabel('temperature', color = 'y')
    # ax3.tick_params('y', colors = 'y')

    # fig.tight_layout()
    # plt.show()

    # fig2= plt.figure()
    # ax4 = fig2.add_subplot(111)
    # # color = plt.cm.viridis()
    # # ax4.scatter(x =np.array(y),y =np.array(hx)  / np.array(y), c=x, cmap = 'viridis')
    # # ax4.plot(np.array(y), np.array(hx)  / np.array(y), linewidth = .5, alpha = .5)
    # # ax4.set_ylabel('hexagonal disorder', color='r')
    # # ax4.set_xlabel('Energy', color='b')
    # # #ax4.set_zlabel('time (epochs)')
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
    print("last Energy is %s"%(y[-1]))
    display_last(best_config, num_lines, canvas, save_best = better, save_name = save_name, draw_lines = True, interact = interact)


# def interact(lattice, num_lines, num_ticks = 4, anneal_epochs = 5000):
#     best_init_config, lattice_type = Best_Random_Config(lattice(num_lines), 
#         lambda: lattice(num_lines),
#         lambda spin, spins_dict: Energy(spin, spins_dict, J = J, H = H),
#         1)

#     global root
#     root = Tk()
#     global canvas
#     canvas = Canvas(root, width = window_size, height= window_size)
#     best_config = sweep_update(best_init_config, 
#         lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift), 
#         num_sweeps_left =1, cur_epochs= 1, total_epochs = 2, show=show, num_ticks = num_ticks)

#     E_avg = Total_Energy(best_init_config, Energy)

#     T = 0 # for now

#     global temp
#     temp = 1000
#     global Eng
#     Eng = E_avg
#     t = tqdm.trange(epochs, desc = "Start")
#     for i in t:
#         best_config =  sweep_update(best_init_config, 
#         lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift),
#         num_sweeps_left = sweeps_per_epoch, cur_epochs = i + 1, total_epochs = epochs, show =show, 
#         former_best_energy = former_best_energy, num_ticks = num_ticks)
#         Eng = Total_Energy(best_config, lambda spin, spins_dict: Energy(spin, spins_dict, J = J))
#         alignment_hex = hex_align_of_hxsqrtr_lattice(best_config)
#         alignment_ng = neighbor_align(best_config)
#         if Eng < former_best_energy:
#             better = True
#         else:
#             better = False
#         #display_during(best_config, num_lines, canvas, save = better, save_name = save_name, draw_lines = True, interact = interact)
#         y.append(Eng)
#         x.append(i)
#         hx.append(alignment_hex)
#         ng.append(alignment_ng)
#         temp = exp_temp(i, epochs)
#         t.set_description("T = %s, E = %s"%("{:0<14f}".format(temp), "{:0<8f}".format(Eng)))
#         data = update_data(data, best_config)
#         T.append(temp)  
#         #print("T = %s, E = %s"%(round(temp, 3), round(Eng, 3)))

#     save_data(data, neighbor_chart, T, y, 'lattice_data/', '%s_%s-pts_%s-epochs'%(lattice_type, len(best_config.keys()), epochs))




#print(Total_Energy(perfect_antif(), lambda spin, spins_dict: Energy(spin, spins_dict, J = J)))
# for tm in tqdm.tqdm(range())1:
#     optimize_and_draw(hexagonal_lattice, num_lines, save_name= "hexagonal_ground_state", epochs = 1000)
#optimize_and_draw(square_lattice, num_lines, save_name= "square_ground_state", epochs = 1000)
# for tm in tqdm.tqdm(range(60)):
#     optimize_and_draw(triangular_lattice, num_lines = 16, save_name= "a_triangular_ground_state_16", epochs = 1000, interact = True)

# optimize_and_draw(snub_square_lattice, num_lines = 10, save_name= "snub_square_ground_state_10", epochs = 6000)
# for tm in tqdm.tqdm(range(40)):
#     optimize_and_draw(hxsqtr_lattice, num_lines = 8, save_name= "hexagon-square-traing_ground_state_8", epochs = 10000, interact = True)
# for tm in tqdm.tqdm(range(1)):
#     optimize_and_draw(hxsqtr_lattice, num_lines = 16, save_name= "hexagon-square-traing_ground_state_16", epochs = 10000, interact = True)
#print(check_all_neighbors_mutual(hxsqtr_lattice()))
kb = 1
sweeps_per_epoch = 1
linger =0
linger_temp = 0.05
temp_exponent = 1


#optimize_and_draw(lambda num_spins: trunc_hex_lattice(spin, num_spins, window_size), num_lines = 4, save_name= "truncated_hexagong_ground_state_4", epochs = 500, num_ticks =60, show = True)

#optimize_and_draw(lambda num_spins: triangular_lattice(spin, num_spins, window_size), num_lines = 8, save_name= "traing_ground_state_4", epochs = 500, num_ticks =60, show = True)

#optimize_and_draw(lambda num_spins: kagome_lattice(spin, num_spins, window_size), num_lines = 6, save_name= "kagome_ground_state_6", epochs = 500, num_ticks =60, show = True)

#optimize_and_draw(lambda num_spins: snbtrhx_lattice(spin, num_spins, window_size), num_lines = 4, save_name= "snbtrhx_ground_state_8", epochs = 500, num_ticks =60, show = True)

#optimize_and_draw(lambda num_spins: hexagonal_lattice(spin, num_spins, window_size), num_lines = 6, save_name= "hexagonal_ground_state_8", epochs = 500, num_ticks =60, show = True)
#optimize_and_draw(lambda num_spins: hxsqtr_lattice(spin, num_spins, window_size), num_lines = 4, save_name= "hexagon-square-traing_ground_state_4", epochs = 100000, num_ticks =60, show = False)

optimize_and_draw(lambda num_spins: square_lattice(spin, num_spins, window_size), num_lines = 12, save_name= "square_ground_state_4", epochs = 1000, num_ticks =20, show = True)
