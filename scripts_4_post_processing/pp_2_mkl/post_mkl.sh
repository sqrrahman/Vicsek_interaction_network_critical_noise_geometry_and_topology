#!/bin/bash

#SBATCH --account=ant_colonies

#SBATCH --partition=normal_q

#SBATCH --nodes=1

#SBATCH --time=0:10:00

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=18






module reset

g++ -O3 -std=c++17 -fopenmp post_mkl.cpp -o post_mkl

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./post_mkl

