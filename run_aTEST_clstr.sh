#!/bin/bash
#SBATCH --job-nam=aTEST
#SBATCH --array 1-15	# Here insert number of subjects
#SBATCH -p cubric-default
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH -o aTEST_%j.out
#SBATCH -e aTEST_%j.err
matlab -nodisplay -nosplash -nodesktop -r "aTEST(${SLURM_ARRAY_TASK_ID});exit;"
