#!/bin/bash
#SBATCH --account=joharav0
#SBATCH --job-name=diag_smm
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joharav@umich.edu

module purge
module load julia
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs Diag

echo "Starting diagnostics at $(date) on $(hostname)"
julia --project -e '
include("diagnostics_runner.jl");
main()
'
echo "Finished at $(date)"
