#!/bin/bash
#SBATCH --job-name=LSSTxSO_cov
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array 3001-4000
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=6:00:00

#SBATCH --account=cosmo
#SBATCH --partition=standard
#SBATCH --qos=user_qos_timeifler
#SBATCH --output=%A.out
#SBATCH --error=%A.err

cd $SLURM_SUBMIT_DIR

module load gsl/2.6
module load python/3.6/3.6.5

./compute_covariances_fourier ${SLURM_ARRAY_TASK_ID} 2 >&/home/u1/xfang/output/job_output_${SLURM_ARRAY_TASK_ID}.log
