from tkinter import *
import numpy as np
import math
from sklearn.preprocessing import normalize 
import copy as copy 


window_size = 512.0
def draw_arrow(canvas, xy0, xy1, tag="changes"):
        canvas.create_line(xy0[0], xy0[1], xy0[0] + xy1[0], xy0[1]+ xy1[1], fill="black", tag=tag)
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
                                                    fill="black", tag=tag)

def norm(d1_array):
    return d1_array / size(d1_array)

def size(d1_array):
    return np.sum(d1_array ** 2) ** .5

def cross_mag(vec1, vec2):
    return (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])

#print(list(range(0, 8)))
    
def string_to_array(string):
    list_chars = []
    for char in string:
        list_chars.append(char)
    list_list_int_chars = [[]]
    for char in string:
        if char in [']', '[', '', ',']:
            pass
        elif char == ' ':
            list_list_int_chars.append([])
        else:
            list_list_int_chars[-1].append(char)
    list_nums = [list_chars_to_num(int_char_list) for int_char_list in list_list_int_chars]
    return list_nums

def list_chars_to_num(list):
    string = ''
    for char in list:
        string += char
    return int(string)

def arr_of_arrs_to_arr_of_arrs_of_arrs(arr_arrs):
    dict_of_hex = {}
    for arr in arr_arrs: 
        if str(arr[:-1]) in dict_of_hex.keys():
            dict_of_hex[str(arr[:-1])].append(arr)
        else:
            dict_of_hex[str(arr[:-1])] = [arr]
    return list(dict_of_hex.values()), dict_of_hex


#print(arr_of_arrs_to_arr_of_arrs_of_arrs([[1, 2, 3], [1, 3, 2], [1, 2, 5], [2, 2, 5], [1, 3, 5]]))