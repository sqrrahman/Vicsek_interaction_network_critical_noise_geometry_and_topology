#!/bin/bash

#SBATCH --account=ant_colonies

#SBATCH --partition=normal_q

#SBATCH --nodes=1

#SBATCH --time=0:10:00

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=9

#SBATCH --mem=1G




module reset

g++ -O3 -std=c++17 -fopenmp post_cc.cpp -o post_cc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./post_cc

