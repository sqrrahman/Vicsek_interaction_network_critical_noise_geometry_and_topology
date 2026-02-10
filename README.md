This repository contains codes and scripts for studying geometry and network topology of Vicsek clusters. The workflow consists of the following stages:

1. Running Vicsek simulations up to steady state,
2. Generating snapshots after reaching steady state,
3. Identifying and analyzing the clusters formed in the snapshots,
4. Post-processing the sampled snapshots,
5. Plotting figures based on the processed data.

Additional scripts are provided for generating animations/movies (scripts_1a_movie) and for plotting the largest clusters (scripts_3a_cluster_shapes).

The codes here are configured for N = 512. User should adjust (a) required number of timesteps to reach steady state, (b) the number of simulations/realizations, (c) the corresponding critical noise values, (d) sampling frequency, (e) number of snapshots/samples etc. as required.

## Contents
Each folder contains required c++ and/or MATLAB files. If a template is provided, use the placeholders to create files from them.
Sample SLURM job scripts and relevant bash command files (commands.txt) are provided to illustrate the typical workflow.

## Requirements
The provided bash scripts are intended for Linux-based HPC systems.
On Windows, the C++ codes can be compiled using a standard compiler (e.g., g++ with OpenMP support) and executed locally (e.g., using WSL, MinGW, or Visual Studio).
On windows, MATLAB scripts can be run interactively.



The overall workflow is divided into the following stages.

### Part 1: Simulation up to steady state
Copy the files in project directory and use the commands provided in commands.txt to generate C++ and job script files and run the jobs.
The noise values should be near critical noise (or the noise value of interest). Adjust noise values until critical noise values are found.
In the project folder, 9 subfolders will be created (named set_1, set_2, ..., set_9) corresponding to 9 density-speed pairs. the subfolders will have the order parameter values and simulation data (position and heading of each particle) at the last timestep.

### Part 1a: Generating movie
Copy/move the files of part_1_simulation_for_movie subfolder in project directory and use the commands provided in commands.txt. Then do the same for part_2_movie_plotting.
The noise values should be the critical noise values (or the noise value of interest).
The movie can be found in project_directory/set_x/movie subfolder.

### Part 2: Sampling (generating snapshots) after steady state
Copy/move the files in project directory and use the commands provided in commands.txt.
Please ensure that the noise values (used in commands) are the corresponding critical noise.
The simulation data at last timestep will be used to do further simulation. The position and heading of each particle will be saved every fixed number of timesteps in 'set_x/sampled' subfolder.

### Part 3: Cluster analysis
Copy/move the files for each folder (such as scripts_3a_mkl) in project directory and use the commands provided in commands.txt. As some of the files will use the previously generated files, run them in order (i.e. scripts_3a, then scripts_3b, ...).
Data from 'sampled' folder will be used for finding the clusters and calculating different parameters. The outputs will be saved in 'sampled' subfolder.

### Part 3a: Plotting cluster shapes
Use the commands provided in commands.txt to run the c++ file and then the matlab file from project_directory/scripts_3a_cluster_shapes folder.
Adjust the matlab code for better visualization.

### Part 4: Post-processing
Use the commands provided in commands.txt to run the jobs from project_directory/script_4_post_processing folder. No need to copy/move files for this one.
It will calculate average over snapshots/samples for each density-speed pair and save them in subfolders (e.g. cmpd_binned in pp_1_cmpd, mkl_binned in pp_2_mkl).

### Part 5: Plotting
Copy the cmpd_binned and similar subfolders and paste them in project_directory\script_5_plotting\post_binned_data\N512 folder.
Run the matlab files in script_5_plotting\plotting_scripts subfolder to generate plots.
