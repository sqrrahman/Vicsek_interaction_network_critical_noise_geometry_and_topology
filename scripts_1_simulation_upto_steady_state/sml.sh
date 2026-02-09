#!/bin/bash

#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=2G

module reset
g++ -O3 -std=c++17 -fopenmp sml.cpp -o sml
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./sml

