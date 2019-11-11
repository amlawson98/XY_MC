from tkinter import *
import numpy as np 
import random
import copy
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex, Normalize
from math import pi as pi
import math
from PIL import ImageGrab
import os
import pickle
import tqdm
from numba import jit
from particle_utils import draw_arrow, size, string_to_array
from lattices import hxsqtr_lattice, triangular_lattice, kagome_lattice, snbtrhx_lattice, hexagonal_lattice, square_lattice, trunc_hex_lattice
from recursive_MC_update import recursive_update

#cmap = plt.get_cmap('gist_rainbow')
window_size = 500.0 
norm = Normalize(vmin=-6, vmax=6)
cmapp = cm.ScalarMappable(cmap='bwr', norm=norm)
cmap = cmapp.get_cmap()
tau = .1
T_0 = 25
end_temp = 0.001
def exp_temp(epochs, epochs_total):
    return T_0* np.exp(-1*(epochs/(-epochs_total/(np.log(end_temp/T_0)*tau))) / tau)



class spin(object):
    def __init__(self, name,  neighbors, pos, direction= random.uniform(0,2*pi), lattice_size = 600):
        self.name = name
        self.direction = direction
        self.neighbors = neighbors
        self.pos = pos

        self.color = 'light blue'
        self.magnitude = lattice_size/15

    def show(self):
        print(self.name, self.direction, self.neighbors, self.pos)

    def radius(self):
        return self.magnitude/6

    def get_corners(self):
        xy0 = np.array(
            [self.pos[0] - self.radius(), self.pos[1] + self.radius()])
        xy1 = np.array(
            [self.pos[0] + self.radius(), self.pos[1] - self.radius()])
        return xy0, xy1

    def mag_moment_xy(self):
        return self.magnitude * np.array([np.cos(self.direction), np.sin(self.direction)])

    def draw(self, canvas, spins_dict, Hamiltonian_particle):
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
        draw_arrow(canvas, self.pos - (self.mag_moment_xy() / 2), self.mag_moment_xy(), tag="changes")
        energy = self.energy(spins_dict, Hamiltonian_particle)

        canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1], tag="changes", fill= rgb2hex(cmapp.to_rgba(energy)[:-1]))

    def draw_lines(self, spins_dict):
        for name in self.neighbors:
            if size(self.pos - spins_dict[name].pos) < (window_size / 2):
                midpoint = self.pos + (spins_dict[name].pos - self.pos) / 2
                canvas.create_line(self.pos[0], self.pos[1], midpoint[0], midpoint[1], fill= 'grey', tag="static")


    def energy(self, others, Hamiltonian_particle):
        return Hamiltonian_particle(self, others, 0)

def display_during(spins_dict, num_lines, canvas, Energy_func, save_name = "dummy_name", time = 10000, draw_lines = True, interact = True, former_best_energy=0):
    for k in spins_dict.keys():
        spins_dict[k].draw(canvas, spins_dict, Energy_func)
        if draw_lines:
            spins_dict[k].draw_lines(spins_dict)
        canvas.pack()
    canvas.update()


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

linger = 0
linger_temp = 0.05
temp_exponent = 1


kb = 1
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
    update_step = 1, i = 1, show = False, former_best_energy = 0, num_ticks =4, draw_lines=True):
    if num_sweeps_left == 0:
        if show:
            canvas.delete("changes")
            display_during(spins_dict, num_lines, canvas, Energy_func, former_best_energy = former_best_energy, draw_lines=draw_lines) 
            #canvas.pack()
            #root.update()
        return spins_dict
    else:
        temp =exp_temp(cur_epochs, total_epochs)
        all_keys = list(spins_dict.keys())
        some_keys = random.choices(all_keys, k=len(all_keys))
        #random.sample(all_keys, math.floor(((num_lines ** 2) * amount_flipped)))
        for spin_key in some_keys:
            cur_spin = spins_dict[spin_key]
            update_spin_statistical(cur_spin, spins_dict, Energy_func, temp, num_ticks)
        return sweep_update(spins_dict, Energy_func, num_sweeps_left - 1, cur_epochs, total_epochs,
         update_step, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks, draw_lines=draw_lines)


J = 1
H = 0
def save_data(data, neighbor_chart, T, E, place, name):
    with open(place + name + '.pkl', 'wb') as f: 
        pickle.dump((data, neighbor_chart, T, E), f, pickle.HIGHEST_PROTOCOL)

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

sweeps_per_epoch = 2
def optimize_and_draw(lattice, num_lines, updater, num_ticks = 4, save_best = True,
 epochs = 5000, interact = True, show = False):
    best_init_config, lattice_type = Best_Random_Config(lattice(num_lines), 
        lambda: lattice(num_lines),
        lambda spin, spins_dict: Energy(spin, spins_dict, J = J, H = H),
        1)
    global root
    root = Tk()
    global canvas
    canvas = Canvas(root, width = window_size, height= window_size)
    best_config= updater(best_init_config, 
        lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift), 
        num_sweeps_left =1, cur_epochs= 1, total_epochs = 2, show=show, draw_lines= True, num_ticks = num_ticks)

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
        best_config= updater(best_init_config, 
        lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = J, H = H, main_theta_shift = theta_shift),
        num_sweeps_left = sweeps_per_epoch, cur_epochs = i + 1, total_epochs = epochs, show =show, num_ticks = num_ticks, draw_lines = False)
        Eng = Total_Energy(best_config, lambda spin, spins_dict: Energy(spin, spins_dict, J = J))
        temp = exp_temp(i, epochs)
        t.set_description("T = %s, E = %s"%("{:0<14f}".format(temp), "{:0<8f}".format(Eng)))
        y.append(Eng)
        x.append(i)
        data = update_data(data, best_config)
        T.append(temp)
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")  
    save_data(data, neighbor_chart, T, y, 'lattice_data/', '%s_%s-pts_%s-J_%s-epochs_%s-ticks_%s-Tmax_%s-Tmin'%(lattice_type, len(best_config.keys()),J, epochs, num_ticks, T_0, end_temp))

def recursive_updater(spins_dict, Energy_func, num_sweeps_left =1, cur_epochs= 1, total_epochs = 2, show=True, draw_lines= True, num_ticks = 48):
    return recursive_update(spins_dict, Energy_func, num_sweeps_left, cur_epochs, total_epochs, exp_temp, kb, split_prob = 0.5, 
    update_step=1, i = num_sweeps_left, show= show, num_ticks =num_ticks, draw_lines= draw_lines)
#optimize_and_draw(lambda num_spins: hxsqtr_lattice(spin, num_spins, window_size), num_lines = 4, save_name= "hexagon-square-traing_ground_state_4", epochs = 500, num_ticks =60, show = True)
#optimize_and_draw(lambda num_spins: square_lattice(spin, num_spins, window_size), num_lines = 8,  epochs = 1000, num_ticks =120, show = False)

#optimize_and_draw(lambda num_spins: triangular_lattice(spin, num_spins, window_size), num_lines = 8, save_name= "traing_ground_state_4", epochs = 50000, num_ticks =60, show = False)

optimize_and_draw(lambda num_spins: kagome_lattice(spin, num_spins, window_size), num_lines = 2, updater = recursive_updater, epochs = 1000, num_ticks =120, show = False)


#optimize_and_draw(lambda num_spins: snbtrhx_lattice(spin, num_spins, window_size), num_lines = 4, save_name= "snbtrhx_ground_state_8", epochs = 50000, num_ticks =60, show = False)

#optimize_and_draw(lambda num_spins: hexagonal_lattice(spin, num_spins, window_size), num_lines = 6, save_name= "hexagonal_ground_state_8", epochs = 500, num_ticks =60, show = True)
