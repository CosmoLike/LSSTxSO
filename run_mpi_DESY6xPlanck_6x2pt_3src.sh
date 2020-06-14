#!/bin/bash
#PBS -S /bin/bash
#PBS -W group_list=cosmo
#PBS -q qualified
### Set the number of nodes,cores and memory that will be used for this job
### select=3 is the node count, ncpus=28 are the cores in each node,
### mem=168gb is memory per node, pcmem=6gb is the memory per core - optional
#PBS -l select=40:ncpus=28:mem=168GB
#PBS -l place=free:shared
#PBS -l cput=5600:00:00
#PBS -l walltime=8:00:00
#PBS -N y3xplanck_mpi
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/


cd $PBS_O_WORKDIR

module load python/2
module load mpich/ge/gcc/64/3.2.1
module load openmpi

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 560 python runDESY6xPlanck_6x2pt_sys_3src.py
date