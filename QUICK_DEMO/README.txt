Monte Carlo Simulation of XY Spins 

NB: To run these programs you will need:

1) Tkinter - This is the animal package I'm usuing. It should come with a normal anaconda distribution. It doesn't work in Jupyter, so make sure to run it on you computer.

2) tdqm - This makes the progress bar, which i think is pretty much essential since it tells you the Temperature and Energy as well as the system anneals. 
Use pip install tqdm. I doubt it'll conflict with any installed packages. Its pretty lightweight. 

A Quick Explainer for what everything is:

	XY_MC.py is the main program here. It centers around a function called optimize_and_draw, that runs the simulated annealing algorithm