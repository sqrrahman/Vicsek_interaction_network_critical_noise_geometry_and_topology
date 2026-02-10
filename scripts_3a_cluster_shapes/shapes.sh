#!/bin/bash

#SBATCH --account=ant_colonies
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --job-name=shapes_analysis
#SBATCH --output=shapes_analysis_%j.out
#SBATCH --error=shapes_analysis_%j.err

module reset

g++ -O3 -std=c++17 -fopenmp shapes.cpp -o shapes

for set_num in {1..9}; do
    echo "Processing set_${set_num}..."
    ./shapes ${set_num}
    echo "Completed set_${set_num}"
done

echo "All cluster analysis complete!"