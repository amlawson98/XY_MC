from tkinter import *

window_size = 400.0

root = Tk()
canvas = Canvas(root, width = window_size, height = window_size)
while True:
	canvas.create_oval(10, 10, 10, 10)
	canvas.pack()