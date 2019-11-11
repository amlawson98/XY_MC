import random
from particle_utils import string_to_array
import matplotlib.pyplot as plt
import numpy as np 
import tqdm
import copy
pi = np.pi
from lattices import kagome_lattice



def square_neigbor_chart(sqroot_num):
    chart = {}
    for i in range(sqroot_num):
        for j in range(sqroot_num):
            neighbors = []
            if i !=0:
                neighbors.append(str([i -1, j]))
            if j !=0:
                neighbors.append(str([i, j-1]))
            if i != (sqroot_num - 1):
                neighbors.append(str([i+1, j]))
            if j != (sqroot_num - 1):
                neighbors.append(str([i, j+1]))
            chart[str([i,j])] = neighbors
    return chart

class marked_vertex(object):
    def __init__(self, neighbors, lattice_size = 10, pos = np.array([0.0, 0.0]), mark = 'none'):
        self.mark = mark
        self.lattice_size = lattice_size
        self.pos = pos
        self.neighbors = neighbors


def mark_all(marked_neighbor_chart, already_marked_keys):
    if len(already_marked_keys) == len(marked_neighbor_chart.keys()):
        return marked_neighbor_chart
    else:
        new_marked_keys = []
        for marked_key in already_marked_keys:
            marked_key_neighbors = marked_neighbor_chart[marked_key].neighbors
            #print("length of marked_key_neighbors is %s"% len(marked_key_neighbors))
            rand_int = random.randint(0,len(marked_key_neighbors) - 1)
            # array_key = string_to_array(marked_key)
            # if rand_int==0:
            #     adj_key = str([array_key[0] + 1, array_key[1]])
            # if rand_int==1:
            #     adj_key = str([array_key[0], array_key[1] + 1])
            # if rand_int==2:
            #     adj_key = str([array_key[0] - 1, array_key[1]])
            # if rand_int==3:
            #     adj_key = str([array_key[0], array_key[1] - 1])
            adj_key = marked_key_neighbors[rand_int]
            if (adj_key in marked_neighbor_chart.keys()):
                if (adj_key not in already_marked_keys + new_marked_keys):
                    marked_neighbor_chart[adj_key].mark = marked_neighbor_chart[marked_key].mark 
                    new_marked_keys.append(adj_key)
        print("size of marked_neighbor_chart is %s, already_marked_keys is %s, new_marked_keys is %s" % (len(marked_neighbor_chart.keys()), len(already_marked_keys), len(new_marked_keys)))
        return mark_all(marked_neighbor_chart, already_marked_keys + new_marked_keys)


def patches(neighbor_chart):
    keys = list(neighbor_chart.keys())
    rand_keys = random.sample(keys, len(keys))
    start1, start2 = rand_keys[0], rand_keys[1]
    marked_vertices = {}
    marked_vertices[start1] = marked_vertex(neighbor_chart[start1], mark='x')
    marked_vertices[start2] = marked_vertex(neighbor_chart[start2], mark='o')
    for key in rand_keys[2:]:
        marked_vertices[key] = marked_vertex(neighbor_chart[key])
    all_marked = mark_all(marked_vertices, [start1, start2])
    x_dict = {}
    o_dict = {}
    for key in all_marked.keys():
        if all_marked[key].mark == 'x':
            x_dict[key] = all_marked[key].neighbors
        if all_marked[key].mark == 'o':
            o_dict[key] = all_marked[key].neighbors
    return x_dict, o_dict

def recursive_probable_split(neighbor_chart, prob):
    if len(neighbor_chart.keys()) == 1:
        return [neighbor_chart]
    elif random.uniform(0, 1) < prob:
        chart1, chart2 = patches(neighbor_chart)
        return recursive_probable_split(chart1, prob)+ recursive_probable_split(chart2, prob)
    else:
        return [neighbor_chart]

kag_latt = kagome_lattice(marked_vertex)
def catalan(n):
    if n==1:
        return 1
    else:
        total= 0
        for i in range(1, n):
            total += catalan(n-i)*catalan(i)
        return total

def get_depth(i, num_sites):
    if catalan(i) > num_sites:
        return i
    else:
        return get_depth(i+1, num_sites) 
def approx_exp_num_patches(prob, num_sites):
    depth = get_depth(1, num_sites)
    total = 0
    for i in range(1, depth + 1):
        total += i*catalan(i)*(prob**(i-1))*((1-prob)**(i))
    return total

def make_histogram(num_shots_per_prob, prob_step):
    probs = np.arange(0, 1, prob_step)
    nums = []
    for prob in tqdm.tqdm(probs):
        total_num = 0
        for shot in range(num_shots_per_prob):
            total_num +=len(recursive_probable_split(square_neigbor_chart(5), prob))
        nums.append(total_num / num_shots_per_prob)
    plt.plot(probs, nums)
    #plt.plot(probs, np.exp(3*probs))
    plt.show()

"""
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
        all_keys = spins_dict.keys()
        some_keys = all_keys
        #random.sample(all_keys, math.floor(((num_lines ** 2) * amount_flipped)))
        for spin_key in some_keys:
            cur_spin = spins_dict[spin_key]
            update_spin_statistical(cur_spin, spins_dict, Energy_func, temp, num_ticks)
        return sweep_update(spins_dict, Energy_func, num_sweeps_left - 1, cur_epochs, total_epochs,
         update_step, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks, draw_lines=draw_lines)
""" 
def energy_of_patch(partial_dict, old_spins_dict, theta_shift, Energy_func):
    total_energy = 0
    num_border_connections = 0
    for border_key in partial_dict.keys():
        total_energy += Energy_func(partial_dict[border_key], old_spins_dict, theta_shift = theta_shift)
        #num_border_connections += len(partial_dict[border_key].neighbors)
    return total_energy, num_border_connections



def update_patch(patch_keys, old_spins_dict, Energy_func, temp, kb, num_ticks):
    all_non_patch_neighbors = []
    all_border_sites = []
    for key in patch_keys:
        all_in_patch = True
        for neighbor in old_spins_dict[key].neighbors:
            if neighbor in patch_keys:
                continue
            else:
                all_in_patch = False
                all_non_patch_neighbors.append(neighbor)
        if not all_in_patch:
            all_border_sites.append(key)
    all_non_patch_neighbors = np.unique(all_non_patch_neighbors)
    all_border_sites = np.unique(all_border_sites)

    new_partial_dict = {}
    for border_site in all_border_sites:
        new_partial_dict[border_site] = copy.copy(old_spins_dict[border_site])
        new_partial_dict[border_site].neighbors = []
        for neighbor in new_partial_dict[border_site].neighbors:
            if neighbor in all_non_patch_neighbors:
                new_partial_dict[border_site].neighbors.append(neighbor)
    possible_thetas = []
    possible_energies = []
    for j in range(num_ticks):
        main_theta_shift = 0
        j_theta_shift = (j / num_ticks)*2*pi
        possible_thetas.append(j_theta_shift + main_theta_shift)
        possible_energies.append(energy_of_patch(new_partial_dict, old_spins_dict, j_theta_shift, Energy_func)[0])

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

    updated_new_partial_dict = {}
    for spin_key in patch_keys:
        updated_new_partial_dict[spin_key] = copy.copy(old_spins_dict[spin_key])
        updated_new_partial_dict[spin_key].direction += new_theta
    return updated_new_partial_dict

def recursive_update(spins_dict, Energy_func, num_sweeps_left, cur_epochs, total_epochs, temp_func, kb, split_prob = 0.5, 
    update_step=1, i = 1, show= False, former_best_energy = 0, num_ticks =48, draw_lines= True):
    #print("temp_func is %s"% temp_func)
    temp = temp_func(cur_epochs, total_epochs)
    all_keys = spins_dict.keys()
    neighbor_chart = {}
    for key in all_keys:
        neighbor_chart[key] = copy.copy(spins_dict[key].neighbors)
    all_split_charts = recursive_probable_split(neighbor_chart, split_prob)
    new_spins_dict = {}
    for split_chart in all_split_charts:
        new_dict = update_patch(split_chart.keys(), spins_dict, Energy_func, temp, kb, num_ticks)
        for key, val in new_dict.items():
            new_spins_dict[key] = val
    return recursive_update(spins_dict, Energy_func, num_sweeps_left - 1, cur_epochs, total_epochs,
        temp_func, kb, former_best_energy = former_best_energy, show = show, num_ticks = num_ticks, draw_lines=draw_lines)       

#make_histogram(100, 0.05)



