#!/bin/bash
#SBATCH --job-name=lsstso6_6x2
#SBATCH --ntasks=1128
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=30:00:00

#SBATCH --account=cosmo
#SBATCH --partition=standard
#SBATCH --qos=user_qos_timeifler
#SBATCH --output=%A.out
#SBATCH --error=%A.err

cd $SLURM_SUBMIT_DIR

module load gsl/2.6
module load python/3.6/3.6.5
module load openmpi3/3.1.4

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpirun -n 1128 python3 xfang_runLSSTxSO_6x2pt_sys_Y6_1sample.py 1128
date