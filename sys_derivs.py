import sys
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from run_cosmolike_mpp import * 
import yaml


def init_lenssample1(param_file):
	params = yaml.load(open(param_file))
	path ="./"
	if "base_dir" in params:
		path = params['base_dir']

	cov_file = path + params['cov_file']
	data_file = path + params['data_file']
	source_nz = path + params['source_nz']
	lens_nz = path + params['lens_nz']
	mask_file = "none"

	ntomo_source = params['ntomo_source']
	ntomo_lens = params['ntomo_lens']
	ntheta = params['tbins']
	theta_min = params['tbounds'][0]
	theta_max = params['tbounds'][1]
	ggl_cut = params['ggl_overlap_cut']
	runmode ="Halofit"
	probes = "3x2pt"
	initcosmo(runmode)
	initsources(source_nz, ntomo_source)
	initlenses(lens_nz, ntomo_lens, Double10(), Double10(),ggl_cut)
	initbins(ntheta, theta_min, theta_max)
	initprobes(probes)
	initdata_real(cov_file, mask_file, data_file)
	(varied_params, 
		cosmo_min, cosmo_fid, cosmo_max, 
		nuisance_min, nuisance_fid, nuisance_max) = parse_priors_and_ranges(params)
	sys_params =[]
	sys_params += ['shear_m_%d'%i for i in xrange(ntomo_source)]
	sys_params += ['source_z_bias_%d'%i for i in xrange(ntomo_source)]
	sys_params += ['lens_z_bias_%d'%i for i in xrange(ntomo_lens)]

	getattr(nuisance_fid, "bias")[0] = 1.7
	getattr(nuisance_fid, "bias")[1] = 1.7
	getattr(nuisance_fid, "bias")[2] = 1.7
	getattr(nuisance_fid, "bias")[3] = 2.0
	getattr(nuisance_fid, "bias")[4] = 2.0
	ndata = 900
	nuisance_fid.print_struct()
	return sys_params, ndata, cosmo_fid, nuisance_fid

def datav_derivs(sys_params, ndata, cosmo_fid, nuisance_fid, outfile, step_width = 1.0):
	npar = len(sys_params)
	derivs = np.zeros((npar,ndata))
	file1 = "FM_datav2"
	nuisance_sigma = InputNuisanceParams().fiducial_sigma()
	print sys_params
	for n,p in enumerate(sys_params):
		pshort = p[:-2]
		i = int(p[-1])
		##Get values on the grid for each param
		np_var = nuisance_fid
		p0 = getattr(nuisance_fid,pshort)[i]
		dp = getattr(nuisance_sigma,pshort)[i]*step_width
		print("evaluting derivative for parameter %s[%d]=%e, delta = %e "%(pshort,i,p0,dp))
		getattr(np_var, pshort)[i]= p0-2.*dp
		#nuisance_fid.print_struct()
		write_cosmolike_datavector(file1,cosmo_fid,np_var)
		dv_mm = np.genfromtxt(file1)[:,1]

		getattr(np_var, pshort)[i]= p0-1.*dp
		write_cosmolike_datavector(file1,cosmo_fid,np_var)
		dv_m = np.genfromtxt(file1)[:,1]
		
		getattr(np_var, pshort)[i]= p0+dp
		write_cosmolike_datavector(file1,cosmo_fid,np_var)
		dv_p = np.genfromtxt(file1)[:,1]
		getattr(np_var, pshort)[i]= p0+2.*dp
		write_cosmolike_datavector(file1,cosmo_fid,np_var)
		dv_pp = np.genfromtxt(file1)[:,1]
		getattr(np_var, pshort)[i]= p0
		derivs[n,:] = (-dv_pp +8.*dv_p -8.*dv_m+dv_mm)/(12.*dp)
	print n
	ind = np.where(np.abs(derivs)  < 1.e-13)
	derivs[ind] = 0.0
	f = open(outfile, "w")
	for i in range(0,ndata):
		for j in range(0,npar):
			f.write("%e "%(derivs[j,i]))
		f.write("\n")
	f.close()
	return derivs

sys_params, ndata, cosmo_fid, nuisance_fid = init_lenssample1("./yaml/Y3_3x2pt_8Mpc_8Mpc_baseline.yaml")
derivs= datav_derivs(sys_params, ndata, cosmo_fid, nuisance_fid,"derivs_lenssample1.txt")

