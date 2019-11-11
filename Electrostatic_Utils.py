from tkinter import *
import numpy as np
import math
from sklearn.preprocessing import normalize 
import copy as copy 

window_size = 512.0
def draw_arrow(canvas, xy0, xy1):
    canvas.create_line(xy0[0], xy0[1], xy0[0] + xy1[0], xy0[1]+ xy1[1], fill="black")
    length = ((xy1[0]) ** 2 + (xy1[1]) ** 2) ** .5
    theta = math.atan2(xy1[1], xy1[0])
    # print("theta = " + str(theta))
    delt_thet = .1
    scale = .9
    xy2 = [scale * length * math.cos(theta + delt_thet), scale * length * math.sin(theta + delt_thet)]
    xy3 = [scale * length * math.cos(theta - delt_thet), scale * length * math.sin(theta - delt_thet)]
    canvas.create_polygon(xy0[0] + xy1[0], xy0[1]+ xy1[1],
                          xy0[0] + xy2[0], xy0[1]+ xy2[1],
                          xy0[0] + xy3[0], xy0[1]+ xy3[1],
                          fill="black")

def norm(d1_array):
  return d1_array / size(d1_array)

def size(d1_array):
  return np.sum(d1_array ** 2) ** .5

def cross_mag(vec1, vec2):
  return (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])

# print(size(np.array([-3, -4])))
# print(np.arctan2(0, -1))
# print(cross_mag(np.array([1, 1]), np.array([-1, 1])))
# x1 = [1, 2, 3]
# x2 = copy.copy(x1)
# x3 = x1
# x3[1] = 5
# x2[1] = 4
# print(x1)
# print(x2)
# print(x3)

