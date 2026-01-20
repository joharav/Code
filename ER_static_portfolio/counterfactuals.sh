#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=counterfactuals
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
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
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_PROJECT=@.

# Headless plotting (GR)
export GKSwstype=nul

echo "Starting counterfactuals.jl at $(date)"
echo "Host: $(hostname)  Threads: ${JULIA_NUM_THREADS}"

julia --threads=${JULIA_NUM_THREADS} -O3 counterfactuals.jl

echo "Job finished at $(date)"