#!/bin/bash
#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

module reset
module load MATLAB

matlab -nodisplay -nosplash -nodesktop -r "shapes_matlab; exit"