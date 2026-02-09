#!/bin/bash

#SBATCH --account=ant_colonies

#SBATCH --partition=normal_q

#SBATCH --nodes=1

#SBATCH --time=0:10:00

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=18



module reset

g++ -O3 -std=c++17 -fopenmp degree_distribution.cpp -o degree_distribution

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./degree_distribution
