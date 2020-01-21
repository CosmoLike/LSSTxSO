import numpy as np
from subprocess import call
import sys

def get_script(i,j):
	script = """
#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q high_pri
#PBS -J %d-%d
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=pack:shared
#PBS -l walltime=7:00:00
#PBS -N LSSTxSO_cov
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u1/xfang/CosmoLike/LSSTxSO/./compute_covariances_fourier $PBS_ARRAY_INDEX >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log
	"""%(i,j)


i_in = int(sys.argv[1])
i_fi = int(sys.argv[2])

script = get_script(i_in, i_fi)

np.savetxt('xfang_cov.txt', script)
call(['qsub', 'xfang_cov.txt'])
print('Job submitted')