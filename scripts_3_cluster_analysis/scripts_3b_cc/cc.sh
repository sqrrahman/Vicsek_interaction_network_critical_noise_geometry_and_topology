#!/bin/bash

#SBATCH --account=ant_colonies

#SBATCH --partition=normal_q

#SBATCH --nodes=1

#SBATCH --time=0:10:00

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=16

#SBATCH --mem=1G




module reset

module load gcc/11.2.0

g++ -O3 -std=c++17 -fopenmp cc.cpp -o cc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./cc

