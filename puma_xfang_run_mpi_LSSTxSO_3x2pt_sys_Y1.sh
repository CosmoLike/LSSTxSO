#!/bin/bash
#SBATCH --job-name=lsstso1_3x2
#SBATCH --ntasks=1152
#SBATCH --ntasks-per-node=96
#SBATCH --nodes=12
#SBATCH --mem-per-cpu=168gb
#SBATCH --time=20:00:00

#SBATCH --account=cosmo
#SBATCH --partition=standard
#SBATCH --qos=user_qos_timeifler
#SBATCH --output=/home/u1/xfang/output/%A.out
#SBATCH --error=/home/u1/xfang/output/%A.err

module load gsl/2/2.1
module load python/2
module load mpich/ge/gcc/64/3.2.1
module load openmpi

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 1152 python xfang_runLSSTxSO_3x2pt_sys_Y1.py
date