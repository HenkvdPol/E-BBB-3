#!/bin/bash
#SBATCH --job-name=ebbb3
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --mail-type=ALL
##SBATCH --mail-user=your.name@institute.nl
##SBATCH --partition=your-partition
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-2

# Submit from the repository folder (where simulate.R is):
#   mkdir -p logs
#   sbatch Additional_files/slurm/sim.slurm
# Every task runs one line of parameters.txt: horizon, R, start seed, use_cluster.

module purge
module load ALICE/default
module load R/4.4.2-gfbf-2024a

cd "$SLURM_SUBMIT_DIR"
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Additional_files/slurm/parameters.txt)

echo "========================================="
echo "Job array ID:  $SLURM_ARRAY_JOB_ID"
echo "Task index:    $SLURM_ARRAY_TASK_ID"
echo "Parameters:    $PARAMS"
echo "========================================="

Rscript simulate.R $PARAMS
