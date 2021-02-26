#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q high_pri
#PBS -J 5001-9045
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=pack:shared
#PBS -l walltime=10:00:00
#PBS -N LSSTxSO_cov
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u1/xfang/CosmoLike/LSSTxSO/./compute_covariances_fourier_1sample $PBS_ARRAY_INDEX 0 >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log




