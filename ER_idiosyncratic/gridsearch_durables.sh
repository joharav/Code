#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --output=gridsearch_durables.log

module load julia  # Load Julia if needed
julia -p 8 gridsearch_durables.jl  # Run Julia script with 8 processes

