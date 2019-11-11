# Monte Carlo Simulation of XY Spins
Simulator and GUI for Magnetism in 2D XY model in Python 3

Alex Lawson 

In the 2-dimensional XY Model each spin is defined by an angle <img src="https://latex.codecogs.com/svg.latex?\Large&space;\theta_i" title="theta_i" /> , and the total energy is given by a dot product with all nearest neighbors term  plus a dot product with an external field term. 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;H=\sum_{<ij>}Jcos{\theta_i-\theta_j}+\sum_{i}Bcos(\theta_i)" title="Hamiltonian" /> 

(Actual energy scales are nonphysical). State of each spin at each iteration chosen from Bolzmann-Gibbs distribution.
<img src="https://latex.codecogs.com/svg.latex?\Large&space;Prob(E_i)=\frac{1}{Z}e^{\frac{-E_i}{k_BT}},Z=\sum_{i}e^{\frac{-E_i}{k_BT}}" title="Gibbs/Bolz" />

__VERY QUICK VERSION__: DOWNLOAD JUST THE FOLDER "QUICK DEMO"- RUN XY_MC_interact.py - Play with controls. 
(Should work with Anaconda distribution default packages)

General Description in 

1) Tkinter - This is the animal package I'm usuing. It should come with a normal anaconda distribution. It doesn't work in Jupyter, so make sure to run it on you computer.

2) tdqm - This makes the progress bar, which i think is pretty much essential since it tells you the Temperature and Energy as well as the system anneals. 
Use pip install tqdm. I doubt it'll conflict with any installed packages. Its pretty lightweight. 

A Quick Explainer for what everything is:

XY_MC.py is the main program here. It centers around a function called optimize_and_draw, that runs the simulated annealing algorithm
