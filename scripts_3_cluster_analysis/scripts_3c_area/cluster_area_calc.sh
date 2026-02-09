#!/bin/bash

#SBATCH --account=ant_colonies

#SBATCH --partition=normal_q

#SBATCH --nodes=1

#SBATCH --time=0:20:00

#SBATCH --ntasks-per-node=1

#SBATCH --cpus-per-task=16

#SBATCH --mem=4G




module reset

module load gcc/11.2.0

g++ -O3 -std=c++17 -fopenmp cluster_area_calc.cpp -o cluster_area_calc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./cluster_area_calc

