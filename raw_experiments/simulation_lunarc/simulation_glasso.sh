#!/bin/bash
#
# job time, change for what your job requires
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --array=1-5
#
# job name
#SBATCH -J data_simulate2
#
# filenames stdout and stderr - customise, include %j
#SBATCH -o process_%A_%a.out
#SBATCH -e process_%A_%a.out

#file output
#SBATCH --mail-user=jonas.wallin@stat.lu.se
#SBATCH --mail-type=END

# write this script to stdout-file - useful for scripting errors
cat $0

# load the modules required for you program - customise for your program
module load GCC/12.3.0  OpenMPI/4.1.5 R/4.4.1

cp -p estimate_simulated.R $NAISS_TMP

export R_LIBS=~/R-packages-4.4.1
# run the program
# customise for your program name and add arguments if required
PART=${SLURM_ARRAY_TASK_ID:-1}
N_PARTS=5
Rscript estimate_simulated.R FALSE "$PART" "$N_PARTS"
