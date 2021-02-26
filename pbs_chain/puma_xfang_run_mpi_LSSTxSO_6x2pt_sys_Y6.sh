#!/bin/bash
#SBATCH --job-name=lsstso6_6x2
#SBATCH --ntasks=1128
#SBATCH --ntasks-per-node=94
#SBATCH --nodes=12
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=20:00:00

#SBATCH --account=cosmo
#SBATCH --partition=standard
#SBATCH --qos=user_qos_timeifler
#SBATCH --output=/home/u1/xfang/output/%A.out
#SBATCH --error=/home/u1/xfang/output/%A.err

module load gnu8/8.3.0
module load gsl/2.6
module load python/3.6/3.6.5
module load mpich/3.3.1
module load openmpi3/3.1.4

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 1128 python xfang_runLSSTxSO_6x2pt_sys_Y6.py
date