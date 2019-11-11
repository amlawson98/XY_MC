**_Monte Carlo Simulation of XY Spins_** 

Alex Lawson 

Simulates Lattice of Spins in 2D XY model, where each spin in defined by an angle theta_i, and the total energy is given by a dot product with all nearest neighbor term  plus a dot product with external field term. 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;H=\sum_{\langle ij \rangle}Jcos(\theta_i - \theta_j) + \sum_{i}H cos(\theta_i)" title="\Large H=\sum_{\langle ij \rangle}Jcos(\theta_i - \theta_j) + \sum_{i}H cos(\theta_i)" /> 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-\sum_0^{10}\beta_i\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

(Actual energy scales are nonphysical). State of each spin at each iteration chosen from Bolzmann-Gibbs distribution.

VERY QUICK VERSION: DOWNLOAD THE FOLDER "QUICK DEMO"- RUN XY_MC_interact.py - Play with controls.

NB: To run these programs you will need:

1) Tkinter - This is the animal package I'm usuing. It should come with a normal anaconda distribution. It doesn't work in Jupyter, so make sure to run it on you computer.

2) tdqm - This makes the progress bar, which i think is pretty much essential since it tells you the Temperature and Energy as well as the system anneals. 
Use pip install tqdm. I doubt it'll conflict with any installed packages. Its pretty lightweight. 

A Quick Explainer for what everything is:

XY_MC.py is the main program here. It centers around a function called optimize_and_draw, that runs the simulated annealing algorithm
