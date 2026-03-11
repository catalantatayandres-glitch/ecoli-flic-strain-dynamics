#!/bin/bash
#SBATCH --job-name=dada2_script
#SBATCH --output=/users/acatala3/ecoli-flic-strain-dynamics/results/03_dada2_testing/logs/dada2_script_%j.out
#SBATCH --error=/users/acatala3/ecoli-flic-strain-dynamics/results/03_dada2_testing/logs/dada2_script_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G

echo "Starting job on $(date)"
echo "Running on node: $HOSTNAME"

# Load modules
module purge
module load miniforge3/24.11.3-2

# Activate your conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate conda_clean

# Move to project directory
cd /users/acatala3/ecoli-flic-strain-dynamics

# Run the pipeline
Rscript scripts/03_dada2_testing/dada2_script_final.R\
  --input_dir /users/acatala3/first_step/pilot_analysis/samples_pilot/trimmed_samples \
  --out_dir /users/acatala3/ecoli-flic-strain-dynamics/results/03_dada2_testing \
  --n_cores 6 \
  --minLen 650 \
  --maxLen 3000

echo "Job finished on $(date)"
