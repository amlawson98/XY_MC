from tkinter import *
import numpy as np 
import random
import copy
from functools import reduce
from matplotlib import pyplot as plt
from math import pi as pi
from math import sin, cos
from particle_utils import draw_arrow
# window_size = 700.0 



# class spin(object):
#     def __init__(self, name, direction, neighbors, pos, magnitude = 60):
#         self.name = name
#         self.direction = direction
#         self.neighbors = neighbors
#         self.pos = pos

#         self.color = 'light blue'
#         self.magnitude = magnitude

#     def show(self):
#         print(self.name, self.direction, self.neighbors, self.pos)

#     def radius(self):
#         return 10

#     def get_corners(self):
#         xy0 = np.array(
#             [self.pos[0] - self.radius(), self.pos[1] + self.radius()])
#         xy1 = np.array(
#             [self.pos[0] + self.radius(), self.pos[1] - self.radius()])
#         return xy0, xy1

#     def mag_moment_xy(self):
#         return self.magnitude * np.array([np.cos(self.direction), np.sin(self.direction)])

#     def draw(self, canvas):
#         xy0, xy1 = self.get_corners()
#         # k = self.name[7]
#         # if k == "0":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
#         #                    outline="black", fill='orange')
#         #     canvas.pack()
#         # if k == "1":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
#         #                    outline="black", fill='green')
#         #     canvas.pack()
#         # if k == "2":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
#         #                    outline="black", fill='red')
#         #     canvas.pack()
#         # if k == "3":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
#         #                    outline="black", fill='light blue')
#         #     canvas.pack()
#         # if k == "4":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
#         #                    outline="black", fill='purple')
#         #     canvas.pack()
#         # if k == "5":
#         #     canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1], outline = "black", fill='pink')
#         #     canvas.pack()
#         draw_arrow(canvas, self.pos - (self.mag_moment_xy() / 2), self.mag_moment_xy())

# # class spin(object):
# #     def __init__(self, name, direction, neighbors, pos):
# #         self.name = name
# #         self.direction = direction
# #         self.neighbors = neighbors
# #         self.pos = pos

# #         self.color = 'light blue'
# #         self.magnitude = 60

# #     def show(self):
# #         print(self.name, self.direction, self.neighbors, self.pos)

# #     def radius(self):
# #         return 10

# #     def get_corners(self):
# #         xy0 = np.array(
# #             [self.pos[0] - self.radius(), self.pos[1] + self.radius()])
# #         xy1 = np.array(
# #             [self.pos[0] + self.radius(), self.pos[1] - self.radius()])
# #         return xy0, xy1

# #     def mag_moment_xy(self):
# #         return self.magnitude * np.array([np.cos(self.direction), np.sin(self.direction)])

# #     def draw(self):
# #         xy0, xy1 = self.get_corners()
# #         canvas.create_oval(xy0[0], xy0[1], xy1[0], xy1[1],
# #                            outline="black", fill=self.color)
# #         draw_arrow(canvas, self.pos - (self.mag_moment_xy() / 2), self.mag_moment_xy())

# #     def energy(self, others, Hamiltonian_particle):
# #         return Hamiltonian_particle(self, others)

def snub_square_lattice(point, number = 4, window_size = 700):
    points = {}
    name = "Snub_Square"
    for i in range(0, number):
        for j in range(0, number):
            for k in range(0, 8):
                if i == 0: 
                    right_neighbor = (number - 1) 
                else: 
                    right_neighbor = (i - 1)
                if i == (number - 1): 
                    left_neighbor = 0 
                else: 
                    left_neighbor = (i + 1)
                if j == 0:
                    up_neighbor = (number - 1)
                else:
                    up_neighbor = j - 1
                if j == (number - 1):
                    down_neighbor = 0
                else:
                    down_neighbor = j + 1
                s = 1 / (1 + (3 ** .5))
                if k == 0:
                    neighbors = [str([i, up_neighbor, 7]), 
                    str([i, up_neighbor, 6]), 
                    str([i, j, 1]), 
                    str([i, j, 2]), 
                    str([i, j, 3])]
                    pos_shift = np.array([(3 ** .5) * s / 2, 0.0])
                elif k == 1:
                    neighbors = [str([right_neighbor, up_neighbor, 7]), 
                    str([i, up_neighbor, 6]),
                    str([right_neighbor, j, 2]), 
                    str([i, j, 0]), 
                    str([i, j, 3])]
                    pos_shift = np.array([((3 ** .5) + 1) * s / 2, 0.0])
                elif k == 2:
                    neighbors = [str([i, up_neighbor, 7]), 
                    str([left_neighbor, j, 1]),
                    str([left_neighbor, j, 5]), 
                    str([i, j, 0]), 
                    str([i, j, 4])]
                    pos_shift = np.array([0.0, s / 2])
                elif k == 3:
                    neighbors = [str([i, j, 0]), 
                    str([i, j, 1]),
                    str([i, j, 4]), 
                    str([i, j, 5]), 
                    str([i, j, 6])]
                    pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, ((3 ** .5) / 2) * s])
                elif k == 4:
                    neighbors = [str([left_neighbor, j, 5]), 
                    str([i, j, 2]),
                    str([i, j, 3]), 
                    str([i, j, 6]), 
                    str([i, j, 7])]
                    pos_shift = np.array([s / 2, (1/2 + ((3 ** .5) / 2)) * s])
                elif k == 5:
                    neighbors = [str([right_neighbor, j, 2]), 
                    str([right_neighbor, j, 4]),
                    str([right_neighbor, j, 7]), 
                    str([i, j, 3]), 
                    str([i, j, 6])]
                    pos_shift = np.array([((1/ 2) + (3 ** .5)) * s, (1/2 + ((3 ** .5) / 2)) * s])
                elif k == 6:
                    neighbors = [str([i, down_neighbor, 0]), 
                    str([i, down_neighbor, 1]),
                    str([i, j, 3]), 
                    str([i, j, 4]), 
                    str([i, j, 5])]
                    pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, (((3 ** .5)/ 2) + 1)* s])
                elif k == 7:
                    neighbors = [
                    str([left_neighbor, down_neighbor, 1]), 
                    str([i, down_neighbor, 2]),
                    str([i, down_neighbor, 0]), 
                    str([left_neighbor, j, 5]), 
                    str([i, j, 4])]
                    pos_shift = np.array([0.0, (.5 + (3 ** .5)) * s])
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([i + 1/2, j + 1/2]) + pos_shift) * (window_size / number))
    return points, name

def perfect_antif(number = 4, window_size = 700):
    name = {}
    for i in range(0, number):
        for j in range(0, number):

            if i == 0: 
                right_neighbor = (number - 1) 
            else: 
                right_neighbor = (i - 1)

            if i == (number - 1): 
                left_neighbor = 0 
            else: 
                left_neighbor = (i + 1)

            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1

            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1

            if ((i + j) % 2) == 0:
                direction = pi / 2
            if ((i + j) % 2) == 1:
                direction = 3 * (pi / 2)
            spins[str([i, j])] = spin(
                name = str([i, j]),
                direction = direction,
                lattice_size = window_size,
                neighbors = [
                str([i, up_neighbor]),
                str([right_neighbor, j]),
                str([left_neighbor, j]),
                str([i, down_neighbor])],
                pos = np.array([i + 1/2, j + 1/2]) * (window_size / number))
    return spins

def pentagon(point, number = 4, window_size = 700):
    points = {}
    name = 'Pentagon'
    for i in range(0, 5):
        position = np.array([1/2, 1/2]) + .1* np.array([np.cos(i*2*np.pi/5), np.sin(i*2*np.pi/5)])
        points[str([i])] = point(
            name = str([i]),
            #direction = random.uniform(0, 2 * pi),
            neighbors = [
            str([(i + 1) % 5]),
            str([(i + 4) % 5])],
            lattice_size = window_size,
            pos = position * (window_size))
    return points, name

def square_lattice(point, number = 4, window_size = 700):
    points = {}
    name = 'Square'
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                right_neighbor = (number - 1) 
            else: 
                right_neighbor = (i - 1)
            if i == (number - 1): 
                left_neighbor = 0 
            else: 
                left_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1
            points[str([i, j])] = point(
                name = str([i, j]),
                #direction = random.uniform(0, 2 * pi),
                neighbors = [
                str([i, up_neighbor]),
                str([right_neighbor, j]),
                str([left_neighbor, j]),
                str([i, down_neighbor])],
                lattice_size = window_size,
                pos = np.array([i + 1/2, j + 1/2]) * (window_size / number))
    return points, name

def triangular_lattice(point, number = 4, window_size = 700):
    points = {}
    name = 'Triangular'
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                right_neighbor = (number - 1) 
            else: 
                right_neighbor = (i - 1)
            if i == (number - 1): 
                left_neighbor = 0 
            else: 
                left_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1

            if j % 2 == 0:
                x_pos = i + .5
                points[str([i, j])] = point(
                    name = str([i, j]),
                    neighbors = [
                    str([i, up_neighbor]),
                    str([right_neighbor, up_neighbor]),
                    str([left_neighbor, j]),
                    str([right_neighbor, j]),
                    str([right_neighbor, down_neighbor]),
                    str([i, down_neighbor])],
                    pos = np.array([x_pos, j]) * (window_size / number))
            else:
                x_pos = i + 1
                points[str([i, j])] = point(
                    name = str([i, j]),
                    neighbors = [
                    str([i, up_neighbor]),
                    str([left_neighbor, up_neighbor]),
                    str([left_neighbor, j]),
                    str([right_neighbor, j]),
                    str([left_neighbor, down_neighbor]),
                    str([i, down_neighbor])],
                    lattice_size = window_size,
                    pos = np.array([x_pos, j]) * (window_size / number))
    return points, name


def kagome_lattice(point, number = 4, window_size = 700):
    points = {}
    name = 'Kagome'
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                left_neighbor = (number - 1) 
            else: 
                left_neighbor = (i - 1)
            if i == (number - 1): 
                right_neighbor = 0 
            else: 
                right_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1
            if j % 2 == 0:
                x_pos = i + .5
                top_right_block = [right_neighbor, up_neighbor]
                top_left_block = [i, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [i, down_neighbor]
                bottom_right_block = [right_neighbor, down_neighbor]
            else:
                x_pos = i 
                top_right_block = [i, up_neighbor]
                top_left_block = [left_neighbor, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [left_neighbor, down_neighbor]
                bottom_right_block = [i, down_neighbor]
            for k in range(0, 3):
                if k == 0:
                    neighbors = [
                    str(right_block + [2]),
                    str(bottom_right_block + [2]),
                    str(bottom_right_block + [1]),
                    str([i, j, 1])]
                    pos_shift = [1/2, 0]
                if k == 1:
                    neighbors = [
                    str(right_block + [2]),
                    str(top_left_block + [0]),
                    str([i, j, 2]),
                    str([i, j, 0])]
                    pos_shift = [1/4, -1/2]
                if k == 2:
                    neighbors = [
                    str([i, j, 1]),
                    str(top_left_block + [0]),
                    str(left_block + [0]),
                    str(left_block + [1])]
                    pos_shift = [-1/4, -1/2]
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    #magnitude = 400 / number,
                    #direction = random.uniform(0, 2 * pi),
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([x_pos, j + 1/2]) + np.array(pos_shift)) * (window_size / number))
    return points, name

def hexagonal_lattice(point, number = 4, window_size = 700):
    points = {}
    name = 'Hexagonal'
    for i in range(0, number):
        for j in range(0, number):
            if (j + i) % 2 == 0:
                neighbors = [str([(i - 1) % number, (j) % number]), 
                str([(i) % number, (j + 1) % number]), 
                str([(i + 1) % number, (j) % number])]
                y_shift = 1/3
            else:
                neighbors = [str([(i - 1) % number, (j) % number]),
                str([(i) % number, (j - 1) % number]),
                str([(i + 1) % number, (j) % number])]
                y_shift = 0
            x_pos = i
            y_pos = j + y_shift
            points[str([i, j])] = point(
                name = str([i, j]),
                neighbors = neighbors,
                lattice_size = window_size,
                pos = np.array([x_pos, y_pos]) * (window_size / number))
    return points, name

def hxsqtr_lattice(point, number = 4, window_size = 700):
    points = {}
    name = "Rhombitrihexagonal"
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                left_neighbor = (number - 1) 
            else: 
                left_neighbor = (i - 1)
            if i == (number - 1): 
                right_neighbor = 0 
            else: 
                right_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1
            if j % 2 == 0:
                x_pos = i + .5
                top_right_block = [right_neighbor, up_neighbor]
                top_left_block = [i, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [i, down_neighbor]
                bottom_right_block = [right_neighbor, down_neighbor]
            else:
                x_pos = i 
                top_right_block = [i, up_neighbor]
                top_left_block = [left_neighbor, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [left_neighbor, down_neighbor]
                bottom_right_block = [i, down_neighbor]
            for k in range(0, 6):
                if k == 0:
                    neighbors = [
                    str(right_block + [2]),
                    str(top_right_block + [4]),
                    str([i, j, 1]),
                    str([i, j, 5])]
                    pos_shift = [1/3, -1/6]
                if k == 1:
                    neighbors = [
                    str(top_right_block + [3]),
                    str(top_left_block + [5]),
                    str([i, j, 2]),
                    str([i, j, 0])]
                    pos_shift = [0, -.5]
                if k == 2:
                    neighbors = [
                    str([i, j, 1]),
                    str(top_left_block + [4]),
                    str(left_block + [0]),
                    str([i, j, 3])]
                    pos_shift = [-1/3, -1/6]
                if k == 3:
                    neighbors = [
                    str([i, j, 2]),
                    str(left_block + [5]),
                    str(bottom_left_block + [1]),
                    str([i, j, 4])]
                    pos_shift = [-1/3, 1/6]
                if k == 4:
                    neighbors = [
                    str([i, j, 5]),
                    str([i, j, 3]),
                    str(bottom_left_block + [0]),
                    str(bottom_right_block + [2])]
                    pos_shift = [0, 0.5]
                if k == 5:
                    neighbors = [
                    str(right_block + [3]),
                    str([i, j, 0]),
                    str([i, j, 4]),
                    str(bottom_right_block + [1])]
                    pos_shift = [1/3, 1/6]
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    #magnitude = 400 / number,
                    #direction = random.uniform(0, 2 * pi),
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([x_pos, j + 1/2]) + np.array(pos_shift)) * (window_size / number))
    return points, name

def snbtrhx_lattice(point, number = 4, window_size = 700):
    points = {}
    name = "Snub-Trihexagonal"
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                left_neighbor = (number - 1) 
            else: 
                left_neighbor = (i - 1)
            if i == (number - 1): 
                right_neighbor = 0 
            else: 
                right_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1
            if j % 2 == 0:
                x_pos = i + .5
                top_right_block = [right_neighbor, up_neighbor]
                top_left_block = [i, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [i, down_neighbor]
                bottom_right_block = [right_neighbor, down_neighbor]
            else:
                x_pos = i 
                top_right_block = [i, up_neighbor]
                top_left_block = [left_neighbor, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [left_neighbor, down_neighbor]
                bottom_right_block = [i, down_neighbor]
            for k in range(0, 6):
                if k == 0:
                    neighbors = [
                    str(right_block + [3]),
                    str(right_block + [2]),
                    str(top_right_block + [4]),
                    str([i, j, 1]),
                    str([i, j, 5])]
                    pos_shift = [1/3, -1/6]
                if k == 1:
                    neighbors = [
                    str(top_right_block + [4]),
                    str(top_right_block + [3]),
                    str(top_left_block + [5]),
                    str([i, j, 2]),
                    str([i, j, 0])]
                    pos_shift = [0, -.5]
                if k == 2:
                    neighbors = [
                    str([i, j, 1]),
                    str(top_left_block + [4]),
                    str(top_left_block + [5]),
                    str(left_block + [0]),
                    str([i, j, 3])]
                    pos_shift = [-1/3, -1/6]
                if k == 3:
                    neighbors = [
                    str([i, j, 2]),
                    str(left_block + [0]),
                    str(left_block + [5]),
                    str(bottom_left_block + [1]),
                    str([i, j, 4])]
                    pos_shift = [-1/3, 1/6]
                if k == 4:
                    neighbors = [
                    str([i, j, 5]),
                    str([i, j, 3]),
                    str(bottom_left_block + [1]),
                    str(bottom_left_block + [0]),
                    str(bottom_right_block + [2])]
                    pos_shift = [0, 0.5]
                if k == 5:
                    neighbors = [
                    str(right_block + [3]),
                    str([i, j, 0]),
                    str([i, j, 4]),
                    str(bottom_right_block + [1]),
                    str(bottom_right_block + [2])]
                    pos_shift = [1/3, 1/6]
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    #magnitude = 400 / number,
                    #direction = random.uniform(0, 2 * pi),
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([x_pos, j + 1/2]) + np.array(pos_shift)) * (window_size / number))
    return points, name

def snub_square_lattice(point, number = 2, window_size = 700):
    points = {}
    name = 'Snub-Square'
    for i in range(0, number):
        for j in range(0, number):
            for k in range(0, 8):
                if i == 0: 
                    right_neighbor = (number - 1) 
                else: 
                    right_neighbor = (i - 1)
                if i == (number - 1): 
                    left_neighbor = 0 
                else: 
                    left_neighbor = (i + 1)
                if j == 0:
                    up_neighbor = (number - 1)
                else:
                    up_neighbor = j - 1
                if j == (number - 1):
                    down_neighbor = 0
                else:
                    down_neighbor = j + 1
                s = 1 / (1 + (3 ** .5))
                if k == 0:
                    neighbors = [str([i, up_neighbor, 7]), 
                    str([i, up_neighbor, 6]), 
                    str([i, j, 1]), 
                    str([i, j, 2]), 
                    str([i, j, 3])]
                    pos_shift = np.array([(3 ** .5) * s / 2, 0.0])
                elif k == 1:
                    neighbors = [str([left_neighbor, up_neighbor, 7]), 
                    str([i, up_neighbor, 6]),
                    str([left_neighbor, j, 2]), 
                    str([i, j, 0]), 
                    str([i, j, 3])]
                    pos_shift = np.array([(((3 ** .5)/ 2)  + 1) * s, 0.0])
                elif k == 2:
                    neighbors = [str([i, up_neighbor, 7]), 
                    str([right_neighbor, j, 1]),
                    str([right_neighbor, j, 5]), 
                    str([i, j, 0]), 
                    str([i, j, 4])]
                    pos_shift = np.array([0.0, s / 2])
                elif k == 3:
                    neighbors = [str([i, j, 0]), 
                    str([i, j, 1]),
                    str([i, j, 4]), 
                    str([i, j, 5]), 
                    str([i, j, 6])]
                    pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, ((3 ** .5) / 2) * s])
                elif k == 4:
                    neighbors = [str([right_neighbor, j, 5]), 
                    str([i, j, 2]),
                    str([i, j, 3]), 
                    str([i, j, 6]), 
                    str([i, j, 7])]
                    pos_shift = np.array([s / 2, (1/2 + ((3 ** .5) / 2)) * s])
                elif k == 5:
                    neighbors = [str([left_neighbor, j, 2]), 
                    str([left_neighbor, j, 4]),
                    str([left_neighbor, j, 7]), 
                    str([i, j, 3]), 
                    str([i, j, 6])]
                    pos_shift = np.array([((1/ 2) + (3 ** .5)) * s, (1/2 + ((3 ** .5) / 2)) * s])
                elif k == 6:
                    neighbors = [str([i, down_neighbor, 0]), 
                    str([i, down_neighbor, 1]),
                    str([i, j, 3]), 
                    str([i, j, 4]), 
                    str([i, j, 5])]
                    pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, (((3 ** .5)/ 2) + 1)* s])
                elif k == 7:
                    neighbors = [
                    str([right_neighbor, down_neighbor, 1]), 
                    str([i, down_neighbor, 2]),
                    str([i, down_neighbor, 0]), 
                    str([right_neighbor, j, 5]), 
                    str([i, j, 4])]
                    pos_shift = np.array([0.0, (.5 + (3 ** .5)) * s])
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    #direction = random.uniform(0, 2 * pi),
                    #magnitude =  400 / number,
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([i, j]) + pos_shift) * (window_size / number))
    return points, name



def trunc_hex_lattice(point, number = 4, window_size = 700):
    points = {}
    name = 'Truncated-Hexagonal'
    for i in range(0, number):
        for j in range(0, number):
            if i == 0: 
                left_neighbor = (number - 1) 
            else: 
                left_neighbor = (i - 1)
            if i == (number - 1): 
                right_neighbor = 0 
            else: 
                right_neighbor = (i + 1)
            if j == 0:
                up_neighbor = (number - 1)
            else:
                up_neighbor = j - 1
            if j == (number - 1):
                down_neighbor = 0
            else:
                down_neighbor = j + 1
            if j % 2 == 0:
                x_pos = i + .5
                top_right_block = [right_neighbor, up_neighbor]
                top_left_block = [i, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [i, down_neighbor]
                bottom_right_block = [right_neighbor, down_neighbor]
            else:
                x_pos = i 
                top_right_block = [i, up_neighbor]
                top_left_block = [left_neighbor, up_neighbor]
                left_block = [left_neighbor, j]
                right_block = [right_neighbor, j]
                bottom_left_block = [left_neighbor, down_neighbor]
                bottom_right_block = [i, down_neighbor]
            for k in range(0, 6):
                if k == 0:
                    neighbors = [
                    str(bottom_right_block+ [3]),
                    str(bottom_right_block+ [4]),
                    str([i, j, 1])]
                    pos_shift = [.5*cos(-1 *pi/12), -.5*sin(-1 *pi/12)]
                if k == 1:
                    neighbors = [
                    str(right_block + [5]),
                    str([i, j, 2]),
                    str([i, j, 0])]
                    pos_shift = [.5*cos(1 *pi/12), -.5*sin(1 *pi/12)]
                if k == 2:
                    neighbors = [
                    str([i, j, 1]),
                    str(right_block + [5]),
                    str([i, j, 3])]
                    pos_shift = [.5*cos(3 *pi/12), -.5*sin(3 *pi/12)]
                if k == 3:
                    neighbors = [
                    str([i, j, 2]),
                    str(top_left_block + [0]),
                    str([i, j, 4])]
                    pos_shift = [.5*cos(5 *pi/12), -.5*sin(5 *pi/12)]
                if k == 4:
                    neighbors = [
                    str([i, j, 5]),
                    str([i, j, 3]),
                    str(top_left_block + [0])]
                    pos_shift = [.5*cos(7 *pi/12), -.5*sin(7 *pi/12)]
                if k == 5:
                    neighbors = [
                    str(left_block + [2]),
                    str([i, j, 4]),
                    str(left_block + [1])]
                    pos_shift = [.5*cos(9 *pi/12), -.5*sin(9 *pi/12)]
                points[str([i, j, k])] = point(
                    name = str([i, j, k]),
                    #magnitude = 400 / number,
                    #direction = random.uniform(0, 2 * pi),
                    neighbors = neighbors,
                    lattice_size = window_size,
                    pos = (np.array([x_pos, j + 1/2]) + np.array(pos_shift)) * (window_size / number))
    return points, name



def serpinksi_lattice(point, number = 4):
    points = {}



# def perfect_snub_square_lattice(number = 4):
#     spins = {}
#     for i in range(0, number):
#         for j in range(0, number):
#             for k in range(0, 8):
#                 if i == 0: 
#                     right_neighbor = (number - 1) 
#                 else: 
#                     right_neighbor = (i - 1)
#                 if i == (number - 1): 
#                     left_neighbor = 0 
#                 else: 
#                     left_neighbor = (i + 1)
#                 if j == 0:
#                     up_neighbor = (number - 1)
#                 else:
#                     up_neighbor = j - 1
#                 if j == (number - 1):
#                     down_neighbor = 0
#                 else:
#                     down_neighbor = j + 1
#                 s = 1 / (1 + (3 ** .5))
#                 if k == 0:
#                     neighbors = [str([i, up_neighbor, 7]), 
#                     str([i, up_neighbor, 6]), 
#                     str([i, j, 1]), 
#                     str([i, j, 2]), 
#                     str([i, j, 3])]
#                     pos_shift = np.array([(3 ** .5) * s / 2, 0.0])
#                     direction = 0
#                 elif k == 1:
#                     neighbors = [str([right_neighbor, up_neighbor, 7]), 
#                     str([i, up_neighbor, 6]),
#                     str([right_neighbor, j, 2]), 
#                     str([i, j, 0]), 
#                     str([i, j, 3])]
#                     pos_shift = np.array([(((3 ** .5)/ 2)  + 1) * s, 0.0])
#                     direction = 0
#                 elif k == 2:
#                     neighbors = [str([i, up_neighbor, 7]), 
#                     str([left_neighbor, j, 1]),
#                     str([left_neighbor, j, 5]), 
#                     str([i, j, 0]), 
#                     str([i, j, 4])]
#                     pos_shift = np.array([0.0, s / 2])
#                     direction = pi
#                 elif k == 3:
#                     neighbors = [str([i, j, 0]), 
#                     str([i, j, 1]),
#                     str([i, j, 4]), 
#                     str([i, j, 5]), 
#                     str([i, j, 6])]
#                     pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, ((3 ** .5) / 2) * s])
#                     direction = pi
#                 elif k == 4:
#                     neighbors = [str([left_neighbor, j, 5]), 
#                     str([i, j, 2]),
#                     str([i, j, 3]), 
#                     str([i, j, 6]), 
#                     str([i, j, 7])]
#                     pos_shift = np.array([s / 2, (1/2 + ((3 ** .5) / 2)) * s])
#                     direction = 0
#                 elif k == 5:
#                     neighbors = [str([right_neighbor, j, 2]), 
#                     str([right_neighbor, j, 4]),
#                     str([right_neighbor, j, 7]), 
#                     str([i, j, 3]), 
#                     str([i, j, 6])]
#                     pos_shift = np.array([((1/ 2) + (3 ** .5)) * s, (1/2 + ((3 ** .5) / 2)) * s])
#                     direction = 0
#                 elif k == 6:
#                     neighbors = [str([i, down_neighbor, 0]), 
#                     str([i, down_neighbor, 1]),
#                     str([i, j, 3]), 
#                     str([i, j, 4]), 
#                     str([i, j, 5])]
#                     pos_shift = np.array([(1/2 + ((3 ** .5) / 2)) * s, (((3 ** .5)/ 2) + 1)* s])
#                     direction = pi
#                 elif k == 7:
#                     neighbors = [
#                     str([left_neighbor, down_neighbor, 1]), 
#                     str([i, down_neighbor, 2]),
#                     str([i, down_neighbor, 0]), 
#                     str([left_neighbor, j, 5]), 
#                     str([i, j, 4])]
#                     pos_shift = np.array([0.0, (.5 + (3 ** .5)) * s])
#                     direction = pi
#                 spins[str([i, j, k])] = spin(
#                     name = str([i, j, k]),
#                     direction = direction,
#                     magnitude =  400 / number,
#                     neighbors = neighbors,
#                     pos = (np.array([i, j]) + pos_shift) * (window_size / number))
#     return spins

def check_all_neighbors_mutual(point_dict):
    Truth = True
    for key in point_dict.keys():
        for neighborkey in point_dict[key].neighbors:
            if key in point_dict[neighborkey].neighbors:
                Truth = Truth and True
            else:
                point_dict[key].show()
                point_dict[neighborkey].show()
                print("+++++++++++++++++++++")
                Truth = Truth and False
    return Truth

num = 6
# for sp in snub_square_lattice(num).keys():
#     snub_square_lattice(num)[sp].show()
#print(check_all_neighbors_mutual(snub_square_lattice(num)))




