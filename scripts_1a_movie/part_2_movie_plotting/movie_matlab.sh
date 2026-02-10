#!/bin/bash

#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=0:30:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

module reset
module load MATLAB

matlab -nodisplay -r "run('movie.m'); exit"

echo "Finished processing movie.m"