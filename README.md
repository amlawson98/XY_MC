# Monte Carlo Simulation of XY Spins
Simulator and GUI for Magnetism with 2D XY model in Python 3

Alex Lawson 

In the 2-dimensional XY Model for spins on a lattice, each spin is defined by an angle <img src="https://latex.codecogs.com/svg.latex?\Large&space;\theta_i" title="theta_i" /> , and the total energy is given by a dot product with all nearest neighbors term  plus a dot product with an external field term. 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;H=\sum_{<ij>}Jcos(\theta_i-\theta_j)+\sum_{i}Bcos(\theta_i)" title="Hamiltonian" /> 

(Actual energy scales are nonphysical). State of each spin at each iteration chosen from Boltzmann-Gibbs distribution.
<img src="https://latex.codecogs.com/svg.latex?\Large&space;Prob(E_i)=\frac{1}{Z}e^{\frac{-E_i}{k_BT}},Z=\sum_{i}e^{\frac{-E_i}{k_BT}}" title="Gibbs/Bolz" />

__VERY QUICK VERSION__: DOWNLOAD JUST THE WHOLE FILE AS A ZIP - __RUN XY_MC_interact.py__ - Play with controls. 
(Should work with Anaconda distribution default packages)


A Quick __Explainer__ for what everything is:

  - General Description and results in Alex_Lawson_Report.pdf

  - XY_MC.py is the main program here. It centers around a function called optimize_and_draw, that runs the simulated annealing algorithm.

  - XY_MC_interact.py pulls up a Tkinter GUI when run, which allows the user to modify lattice type, temperature and Hamiltonian at will. 
  
  - The rest is all different versions/subprojects, function libraries to import, and diagrams/pictures for presentations and reports.

![Alt text](Notable_XY_Sim_Images/XY_MC_interact_screenshot.PNG?raw=true "Screenshot for XY_MC_interact.py")

__Dependencies__
1) Python 3, with packages from default Anaconda Distribution
    - Using packages: tkinter (For GUI), numpy, random, copy, matplotlib, math, time, numba, sklearn, PIL

2) tdqm - This makes the progress bar. This essential for anything **except XY_MC_interact.py**, since it tells you the Temperature and Energy as well as the system anneals. 
    - Use pip install tqdm. I doubt it'll conflict with any installed packages. Its pretty lightweight. 
