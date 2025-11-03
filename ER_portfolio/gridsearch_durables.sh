#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=smm
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --requeue
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joharav@umich.edu

set -euo pipefail

module purge
module load julia

mkdir -p logs

# Map BLAS/LAPACK to 1 thread to avoid oversubscription with Julia threads
export JULIA_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_PROJECT=@.

# Headless plotting (GR)
export GKSwstype=nul

echo "Starting gridsearch_durables.jl at $(date)"
echo "Host: $(hostname)  Threads: ${JULIA_NUM_THREADS}"

julia --threads=${JULIA_NUM_THREADS} -O3 gridsearch_durables.jl

echo "Job finished at $(date)"
