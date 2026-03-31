#!/bin/bash
#
# job time, change for what your job requires
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-7
#
# job name
#SBATCH -J data_simulate_big
#
# filenames stdout and stderr - customise, include %j
#SBATCH -o process_big_%A_%a.out
#SBATCH -e process_big_%A_%a.out

#file output
#SBATCH --mail-user=jonas.wallin@stat.lu.se
#SBATCH --mail-type=END

module load GCC/12.3.0  OpenMPI/4.1.5 R/4.4.1

export R_LIBS=~/R-packages-4.4.1

PART=${SLURM_ARRAY_TASK_ID:-1}
N_PARTS=7
# Set N_LIST="all" or "200,300,5000" to change the n values.
N_LIST=${N_LIST:-5000}
RESULTS_DIR=${RESULTS_DIR:-results}

mkdir -p "$RESULTS_DIR"

echo "Starting big run: part ${PART} of ${N_PARTS}, n_list=${N_LIST}, results_dir=${RESULTS_DIR}"
OUTPUT_DIR="$RESULTS_DIR" Rscript estimate_simulated_big.R TRUE "$PART" "$N_PARTS" "$N_LIST"
echo "Finished big run: part ${PART} of ${N_PARTS}"
