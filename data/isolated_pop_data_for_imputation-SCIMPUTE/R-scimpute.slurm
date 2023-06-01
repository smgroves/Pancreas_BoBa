#!/bin/bash
#SBATCH --mail-user=christian.t.meyer@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=128:00:00
#SBATCH --mem=16G
#SBATCH --array=1-12
#SBATCH --output=Job_%A_%a.out


echo "Time Started:"
date
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

#Load modules
module load GCC OpenMPI R R-bundle-Bioconductor

Rscript --no-save scImpute-pancreas_slurm.r $SLURM_ARRAY_TASK_ID

#Print finish
date
