#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=evalfun
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=96:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joharav@umich.edu

module purge
module load julia

mkdir -p logs
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Starting evalfun.jl at $(date) on $(hostname)"
echo "CPUs per task: $SLURM_CPUS_PER_TASK | JULIA_NUM_THREADS: $JULIA_NUM_THREADS"

julia --project=. evalfun.jl

echo "Finished evalfun.jl at $(date)"
