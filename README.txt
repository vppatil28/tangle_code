Code for simulating tangling and untangling worms using the Kirchhoff model for elastic fibers.

In this repository:

worm_simulator.m - code for simulating worms tangling and untangling. Requires the following to be in the same directory: 
qualitative12.m
input_elastic_parameters.txt
input_active_head_force_parameters.txt
initial_conditions folder 
input_initial_condition.txt
results folder 

worm_simulator_fast.m - code for simulating worms tangling and untangling at high speed. Contact handling is slightly modified relative to worm_simulator.m. Requires the following to be in the same directory: 
qualitative12.m
input_elastic_parameters_fast.txt
input_active_head_force_parameters_fast.txt
initial_conditions folder
input_initial_condition.txt
results folder 

qualitative12.m - contains colormap used in the worm simulation code

input_elastic_parameters.txt - plaintext file containing elastic input parameters and instructions on how to modify them. This input file is read by worm_simulator.m 

input_elastic_parameters_fast.txt - plaintext file containing elastic input parameters for worm_simulator_fast.m and instructions on how to modify them. This input file is read by worm_simulator_fast.m 

input_active_head_force_parameters.txt - plaintext file containing input parameters controlling the active head force driving the worm and instructions on how to modify them. This input file is read by worm_simulator.m 

input_active_head_force_parameters_fast.txt - plaintext file containing input parameters controlling the active head force driving the worm and instructions on how to modify them. This input file is read by worm_simulator_fast.m 


initial_conditions - folder containing .mat files which store the initial worm positions, and the size of the domain in which they are contained.

input_initial_condition.txt - plaintext file specifying which initial condition will be simulated, and instructions on how to choose the initial condition.

results - folder which will contain output files from running the worm simulation code.


Instructions for using worm_simulator.m and worm_simulator_fast.m

Some choices of initial configurations are given in .mat files in the initial_conditions folder. Each .mat file contains N worms, stored in an array W0 (size 3N x n, where n is the number of discrete points along each worm). W0 represents the spatial positions of the N worms.

A choice of initial condition, material parameters and active head force parameters can be specified by editing plaintext input files (as described above).

Running the worm simulation scripts will simulate worm motion, with the given initial condition, material parameters and active head force parameters. worm_simulator.m and worm_simulator_fast.m both run for 40000 time steps. These parameters can also be edited within the scripts if desired. All quantities are given in simulation units, L (length), M (mass), T (time). To convert to physical units, choose explicit values for L, M, T. 

Output will appear in the form of a video and a .mat file in the results folder. Output data consist of W and TwDispW (along with other auxiliary data described within the code): 
W is a 3N x n x 201 array, representing the spatial position of each worm at 200 evenly spaced timepoints across the simulation. 
TwDispW is an N x (n-2) x 201 array for each worm in the .mat file, representing the twist along each worm at equally spaced timepoints. 
For both position and twist arrays, the last dimension is the time dimension.





