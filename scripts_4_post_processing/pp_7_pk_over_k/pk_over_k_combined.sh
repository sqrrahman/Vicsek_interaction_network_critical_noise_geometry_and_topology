#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:10:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=1G



module reset
g++ -O3 -std=c++17 -fopenmp pk_over_k.cpp -o pk_over_k
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./pk_over_k


module reset
g++ -O3 -std=c++17 -fopenmp pk_combined.cpp -o pk_combined
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./pk_combined

