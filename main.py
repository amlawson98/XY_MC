from tkinter import *
import numpy as np 
import random
import copy
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex, Normalize
from math import pi as pi
import time
import math
from numba import jit
from particle_utils import draw_arrow, size, string_to_array
from lattices import hxsqtr_lattice, triangular_lattice, kagome_lattice, snbtrhx_lattice
from lattices import hexagonal_lattice, square_lattice, trunc_hex_lattice, snub_square_lattice, pentagon


norm = Normalize(vmin=-6, vmax=6)
cmapp = cm.ScalarMappable(cmap='bwr', norm=norm)
cmap = cmapp.get_cmap()


class spin(object):
    def __init__(self, name, neighbors, pos, direction= random.uniform(0,2*pi), current_energy = 0, lattice_size = 600):
        self.name = name
        self.direction = direction
        self.neighbors = neighbors
        self.pos = pos
        self.cur_energy = current_energy
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

    def draw_2(self, canvas, spins_dict):#, Hamiltonian_particle):
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

        #energy = self.energy(spins_dict, Hamiltonian_particle)

        canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1], tag="changes", fill= rgb2hex(cmapp.to_rgba(self.cur_energy)[:-1]))


    def draw_lines(self, spins_dict):
        for name in self.neighbors:
            if size(self.pos - spins_dict[name].pos) < (window_size / 2):
                midpoint = self.pos + (spins_dict[name].pos - self.pos) / 2
                canvas.create_line(self.pos[0], self.pos[1], midpoint[0], midpoint[1], fill= 'grey', tag="user_changes")


    def energy(self, others, Hamiltonian_particle):
        return Hamiltonian_particle(self, others, 0)

def display_during(spins_dict, canvas, Energy_func, save_name = "dummy_name", time = 10000, draw_lines = True, interact = True, former_best_energy=0):
    if draw_lines:
        canvas.delete("user_changes")
    for k in spins_dict.keys():
        spins_dict[k].draw(canvas, spins_dict, Energy_func)
        if draw_lines:
            spins_dict[k].draw_lines(spins_dict)
        canvas.pack()
    canvas.update()


def display_during_2(spins_dict, canvas, save_name = "dummy_name", time = 10000, draw_lines = True, interact = True, former_best_energy=0):
    if draw_lines:
        canvas.delete("user_changes")
    for k in spins_dict.keys():
        spins_dict[k].draw_2(canvas, spins_dict)
        if draw_lines:
            spins_dict[k].draw_lines(spins_dict)
    canvas.update()


def Energy(main_spin, spins_dict, J = 1, H = 0, main_theta_shift = 0):
    directions = list(map(lambda name: spins_dict[name].direction, main_spin.neighbors))
    main_theta = main_spin.direction + main_theta_shift
    return Basic_Energy(main_theta, directions, J, H, main_theta_shift)

@jit
def Basic_Energy(main_theta, other_thetas, J, H, main_theta_shift):
    E = H * np.cos(main_theta + main_theta_shift)
    for theta in other_thetas:
        E += (J * np.cos(main_theta +main_theta_shift - theta))
    return E

def Total_Energy(spins_dict):
    E_tot = 0
    for spin_key in spins_dict.keys():
        E_tot += spins_dict[spin_key].cur_energy
    return E_tot

def Best_Random_Config(spins_dict, generate_spins, num_tries_left):
    if num_tries_left == 0:
        return spins_dict
    else:
        new_spins = generate_spins()
        #print(Total_Energy(spins_dict, Energy_func))
        #print(Total_Energy(new_spins, Energy_func))
        if Total_Energy(spins_dict[0]) <= Total_Energy(new_spins[0]):
            best_spins = copy.copy(spins_dict)
        else:
            best_spins = copy.copy(new_spins)
        #print(Total_Energy(best_spins, Energy_func))
        return Best_Random_Config(best_spins, generate_spins, num_tries_left - 1)


kb = 1
@jit
def scrub_infs(relative_probabilities, possible_energies, temp, E_shift = 0):
    if np.any(np.isinf(relative_probabilities)):
        E_increment = 0.1
        new_relative_probabilities = np.exp((-np.array(possible_energies) - E_shift - 0.1)/(temp * kb))
        return scrub_infs(new_relative_probabilities, possible_energies, temp, E_shift + 0.1) 
    else:
        #print("E_shift = %s"%E_shift)
        return relative_probabilities
    # return relative_probabilities

@jit
def update_spin_fast(main_theta, other_thetas, basic_energy_func, J, H, temp, num_ticks):#, prob_num):
    possible_thetas = []
    possible_energies = []
    for j in range(num_ticks):
        j_theta_shift = (j / num_ticks)*2*pi
        possible_thetas.append(j_theta_shift + main_theta)
        possible_energies.append(basic_energy_func(main_theta, other_thetas, J, H, j_theta_shift))
    relative_probabilities_old = np.exp(-np.array(possible_energies)/(temp * kb))
    relative_probabilities = scrub_infs(relative_probabilities_old, possible_energies, temp, E_shift =0)
    #print(relative_probabilities)
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
    new_energy = possible_energies[index]
    return new_theta, new_energy

def update_spin_statistical_2(spin, spin_dict, basic_Energy_func, J, H, temp, num_ticks):
    main_theta = spin.direction
    other_thetas = list(map(lambda name: spin_dict[name].direction, spin.neighbors))
    new_theta, new_energy = update_spin_fast(main_theta, other_thetas, basic_Energy_func, J, H, temp, num_ticks)#, random.uniform(0,1))
    spin.direction = new_theta
    spin.cur_energy = new_energy

def update_spin_statistical(spin, spin_dict, Energy_func, temp, num_ticks):
    main_theta = spin.direction
    possible_thetas = []
    possible_energies = []
    for j in range(num_ticks):
        j_theta_shift = (j / num_ticks)*2*pi
        possible_thetas.append(j_theta_shift + main_theta)
        possible_energies.append(Energy_func(spin, spin_dict, j_theta_shift))
    relative_probabilities = np.exp(-np.array(possible_energies)/(temp * kb))
    Z = np.sum(relative_probabilities)
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

def sweep_update_2(spins_dict, basic_energy_func, num_sweeps_left, J, H, temp,
    update_step = 1, i = 1, show = False, former_best_energy = 0, num_ticks =4, draw_lines=True):
    if num_sweeps_left == 0:
        if show:
            canvas.delete("changes")
            display_during_2(spins_dict, canvas, former_best_energy = former_best_energy, draw_lines=draw_lines) 
        return spins_dict
    else:
        all_keys = list(spins_dict.keys())
        np.random.shuffle(all_keys)
        for spin_key in all_keys:
            cur_spin = spins_dict[spin_key]
            update_spin_statistical_2(cur_spin, spins_dict, basic_energy_func, J, H, temp, num_ticks)
        return sweep_update_2(spins_dict, basic_energy_func, num_sweeps_left - 1, J, H, temp,
         update_step, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks, draw_lines=draw_lines)

def sweep_update(spins_dict, Energy_func, num_sweeps_left, temp,
    update_step = 1, i = 1, show = False, former_best_energy = 0, num_ticks =4, draw_lines=True):
    if num_sweeps_left == 0:
        if show:
            canvas.delete("changes")
            display_during(spins_dict, canvas, Energy_func, former_best_energy = former_best_energy, draw_lines=draw_lines) 
        return spins_dict
    else:
        all_keys = list(spins_dict.keys())
        some_keys = all_keys
        for spin_key in some_keys:
            cur_spin = spins_dict[spin_key]
            update_spin_statistical(cur_spin, spins_dict, Energy_func, temp, num_ticks)
        return sweep_update(spins_dict, Energy_func, num_sweeps_left - 1, temp,
         update_step, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks, draw_lines=draw_lines)

####################
# BEN THIS IS WHERE TKINTER STUFF HAPPENS
####################

window_width = 800 
window_height = 600
window_size = window_height
side_space= 200
title_width = 60
title_height = 5
root = Tk()
f0 = Frame(root, width = window_width, height=window_height + title_height)
f0.pack(side=TOP, anchor = N)
f1 = Frame(f0, width = window_width - side_space, height = window_height)
f1.pack(side=RIGHT, anchor=W)
f2 = Frame(f0, width =side_space, height = window_height)
f2.pack(side=LEFT, anchor=E)


title = Text(f1, width= title_width, height = title_height)
title.tag_configure('Bold', font = ('Times', 14, 'bold'), justify =CENTER)
title.tag_configure('Big', font = ('Times', 14), justify =CENTER)
title.tag_configure('Small', font = ('Times', 11), justify =CENTER)
title.insert(END, "Monte Carlo Simulation of 2D XY Spin Model \n", 'Bold')
title.insert(END, "H = J\u2211_\u27E8ij\u27E9 cos(\u03B8_i - \u03B8_j) + B\u2211_i cos(\u03B8_i) \n", 'Big')
title.insert(END, "Energy per site given by J * dot product with each nearest neighbor" +
    " plus dot \n product with B field. Periodic boundary conditions on all lattices", 'Small')
title.pack(side = TOP, anchor = N)

canvas = Canvas(f1, width = window_width - side_space, height= window_height - title_height)
lattice_width = window_height - 50

global Avg_E
Avg_E = 0
Energy_per_site_Counter = Text(f2, width = 20, height= 2)
Energy_per_site_Counter.tag_configure('Small', font = ('Times', 11), justify =CENTER)
def update_energy_counter():
    Energy_per_site_Counter.delete(1.0, END)
    Energy_per_site_Counter.insert(END, 'AVG Energy per Site \n = %s'%Avg_E, 'Small')
    Energy_per_site_Counter.pack()
update_energy_counter()

def set_with_temp_scale(scale_val):
    global scale_temp
    scale_temp = np.exp(float(scale_val))

temp_scale = Scale(f2, from_=5, to=-10, length=300, resolution= 0.001, orient= VERTICAL, showvalue=1, command= set_with_temp_scale, label = "log(Temp)")
temp_scale.place()#rely = 30, relx = 80)
temp_scale.pack()#fill=BOTH, side=LEFT)

def set_with_J_scale(scale_val):
    global scale_J
    scale_J = float(scale_val)

J_scale = Scale(f2, from_=1, to=-1,resolution= 1, orient= HORIZONTAL, showvalue=1, command= set_with_J_scale, label = "J (Interaction Energy)")
J_scale.place()#rely = 30, relx = 80)
J_scale.pack(fill=BOTH)

def set_with_H_scale(scale_val):
    global scale_H
    scale_H = float(scale_val)

H_scale = Scale(f2, from_=5, to=0,resolution= 0.01, orient= HORIZONTAL, showvalue=1, command= set_with_H_scale, label = "B (External Field Energy)")
H_scale.place()
H_scale.pack(fill=BOTH)

running = False
txt_var = StringVar()
txt_var.set("Play")

button_frame= Frame(f2)
button_frame.pack()
def pause_play():
    global running
    if running:
        txt_var.set("Play")
        pause_button.configure(bg ="light green")
    if not running:
        txt_var.set("Pause")  
        pause_button.configure(bg = "coral")
    running = not running

pause_button = Button(button_frame, textvariable=txt_var, command= pause_play, bg="light green")
pause_button.pack(side=LEFT)


make_a_step = False
def make_step():
    global make_a_step
    make_a_step = True
step_button = Button(button_frame, text="Step", command=make_step)
step_button.pack(side=LEFT)

quit_button = Button(button_frame, text="Quit", command = root.destroy)
quit_button.pack(side=LEFT)

annealing = BooleanVar()
annealing_button = Checkbutton(f2, text="Annealing", variable=annealing)
annealing_button.pack()
temp=0
scale_J=0
scale_H=0

def update_temp():
    global temp
    decrease_rate = -0.02
    temp = temp*(1 + decrease_rate)

name_string = StringVar()
lattices_dict = {'Square': square_lattice,'Snub-Square': snub_square_lattice, 'Triangular': triangular_lattice,
    'Hexagonal': hexagonal_lattice, 'Kagome': kagome_lattice, 
    'Rhombitrihexagonal': hxsqtr_lattice, 'Snub-Trihexagonal': snbtrhx_lattice,
    'Truncated-Hexagonal': trunc_hex_lattice, 'Pentagon': pentagon}
default_number_dict = {'Square': 10,'Snub-Square': 3, 'Triangular': 10,
    'Hexagonal': 10, 'Kagome': 6, 
    'Rhombitrihexagonal': 6, 'Snub-Trihexagonal': 6,
    'Truncated-Hexagonal': 6, 'Pentagon': 1}
lattice_choices = ["Lattice Type: \n" + key for key in lattices_dict.keys()]


latticeMenu = OptionMenu(f2, name_string, *lattice_choices)# text="Lattice")

# on change dropdown value
def change_dropdown(*args):
    global lattice_type
    string = name_string.get()
    lattice_type = string[15:]  

# link function to change dropdown
name_string.trace('w', change_dropdown)
latticeMenu.pack()


def optimize_and_draw(lattice, num_lines, num_ticks = 4, interact = True, show = False, sweeps_per_epoch=1, Total_Energy_calc = False):
    current_lattice = lattice
    lattice_size_num = num_lines
    global lattice_type
    global temp
    global scale_J
    global scale_H
    global make_a_step
    global canvas
    global best_config
    global Avg_E
    best_config, current_lattice_name = Best_Random_Config(current_lattice(lattice_size_num), 
    lambda: lattice(num_lines),
    1)
    num_sites = len(list(best_config.keys()))
    lattice_type = current_lattice_name
    name_string.set("Lattice Type: \n" + lattice_type)
    latticeMenu.pack()
    scale_J, scale_H = 0, 0 
    best_config= sweep_update(best_config, 
        lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = 0, H = 0, main_theta_shift = theta_shift), 
        num_sweeps_left=1, temp=10, show=show, draw_lines= True, num_ticks = num_ticks)
    while True:
        if annealing.get() and running:
            update_temp()
            temp_scale.set(np.log(temp))
        if lattice_type != current_lattice_name:
            current_lattice = lambda num_spins: lattices_dict[lattice_type](spin, num_spins, lattice_width)
            lattice_size_num = default_number_dict[lattice_type]
            best_config, current_lattice_name = Best_Random_Config(current_lattice(lattice_size_num), 
            lambda: current_lattice(lattice_size_num),
            1)
            num_sites = len(list(best_config.keys()))
            best_config= sweep_update_2(best_config, Basic_Energy,
                num_sweeps_left = 1, J =scale_J, H =scale_H, temp=10, show=show, draw_lines= True, num_ticks = num_ticks)
        if running:
            temp = np.exp(temp_scale.get())
            scale_J = J_scale.get()
            scale_H = H_scale.get()
            #anon_energy_func = lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = scale_J, H = scale_H, main_theta_shift = theta_shift)
            best_config= sweep_update_2(best_config, Basic_Energy,
            num_sweeps_left = sweeps_per_epoch, J =scale_J, H =scale_H, temp=temp, show =show, num_ticks = num_ticks, draw_lines = False)
            if Total_Energy_calc:
                Eng = Total_Energy(best_config)
                Avg_E = Eng/ num_sites
                update_energy_counter()
        if make_a_step:
            temp = np.exp(temp_scale.get())
            #anon_energy_func = lambda spin, spins_dict, theta_shift: Energy(spin, spins_dict, J = scale_J, H = scale_H, main_theta_shift = theta_shift)
            best_config= sweep_update_2(best_config, Basic_Energy,
            num_sweeps_left = sweeps_per_epoch,  J =scale_J, H =scale_H, temp=temp, show =show, num_ticks = num_ticks, draw_lines = False)
            if Total_Energy_calc:
                Eng = Total_Energy(best_config)
                Avg_E = Eng/ num_sites
                update_energy_counter()
            make_a_step= False
        canvas.update()
        f1.update()
    root.mainloop()

optimize_and_draw(lambda num_spins: square_lattice(spin, num_spins, lattice_width), num_lines = 10, num_ticks =300, show = True, Total_Energy_calc = True)
