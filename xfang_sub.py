import numpy as np
from subprocess import call
import sys

def get_script(i,j, scenario):
	script = """
#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q high_pri
#PBS -J %d-%d
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=pack:shared
#PBS -l walltime=5:00:00
#PBS -N LSSTxSO_cov
#PBS -e /home/u1/xfang/output/
#PBS -o /home/u1/xfang/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u1/xfang/CosmoLike/LSSTxSO/./compute_covariances_fourier $PBS_ARRAY_INDEX %d >&/home/u1/xfang/output/job_output_$PBS_ARRAY_INDEX.log
	"""%(i,j, scenario)
	return script


i_in = np.int(sys.argv[1])
i_fi = np.int(sys.argv[2])
scenario = np.int(sys.argv[3])

print(i_in,i_fi)
print('scenario: %d'%(scenario))
script = get_script(i_in, i_fi, scenario)
print(script)

f = open("xfang_cov.txt","w")
f.write(script)
f.close()
call(['qsub', 'xfang_cov.txt'])
print('Job submitted')