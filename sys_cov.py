import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import yaml

sigma_m = np.array([0.005,0.005,0.005,0.005])
sigma_source_pz = np.array([0.005,0.005,0.005,0.005])
sigma_lens_pz = np.array([0.004,0.003,0.003,0.005,0.011])

#modify these to account for correlations across z-bins
prior_shear_m = np.diag(sigma_m**2)
prior_source_dz = np.diag(sigma_source_pz**2)
prior_lens_dz = np.diag(sigma_lens_pz**2)


def build_sys_cov(derivs, prior):
	nsys = np.shape(derivs)[0]
	ndata = np.shape(derivs)[1]
	#print np.shape(derivs), np.shape(prior)
	if (nsys != np.shape(prior)[0]):
		print "build_sys_cov: dimension of derivatives and prior don't match"
		exit
	cov = np.zeros((ndata,ndata))
	for i in range(nsys):
		for j in range(nsys):
			cov +=prior[i,j]*np.tensordot(derivs[i,:],derivs[j,:],axes = 0)
	return cov
def read_cov(covfile):
	data = np.genfromtxt(covfile)
	cov = np.zeros((int(np.max(data[:,0])+1),int(np.max(data[:,0])+1)))
	for i in range(0,data.shape[0]):
  		cov[int(data[i,0]),int(data[i,1])] =data[i,2]
  		cov[int(data[i,1]),int(data[i,0])] =data[i,2]
  	return cov

covfile = "cov_Y3/cov_cosmolike_y3_mcal_DESY1_bf_lens1.txt"
derivs_file = "derivs_lenssample1.txt"
outfile = covfile[:-4]+"_sys.txt"
cov = read_cov(covfile)
derivs = np.transpose(np.genfromtxt(derivs_file))

derivs_m = derivs[0:4,:]
derivs_source_dz = derivs[4:8,:]
derivs_lens_dz = derivs[8:13,:]
ndata = np.shape(derivs)[1]

cov_sys = build_sys_cov(derivs_m, prior_shear_m)
cov_sys += build_sys_cov(derivs_source_dz, prior_source_dz)
cov_sys += build_sys_cov(derivs_lens_dz, prior_lens_dz)

cov += cov_sys
f = open(outfile, "w")
for i in range(0,ndata):
	for j in range(0,ndata):
		f.write("%d %d %e\n"%(i,j,cov[j,i]))
f.close()

plt.subplot(1, 1, 1)
ax = plt.gca()
im = ax.imshow(np.log(cov_sys), interpolation="nearest",origin='lower')
plt.colorbar(im)
plt.show()

