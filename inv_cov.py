#!/usr/bin/python
import sys
import math, numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy import linalg as LA
import numpy as np


infile =['/users/timeifler/Dropbox/cosmolike_store/LSSTxSO/cov/cov_LSSTY10xCMBS4']

#infile =['/users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/cov/cov_WFIRST_Ncl25_4clusterbins_nrichmin25_source_Dec17']
data= ['datav/6x2pt_LSST_Y10_area=1.800000e+04_dmo']
outname=['cov_LSSTY10xCMBS4']

# the numbers below can be computed knowing the data vector settings, e.g. 10 tomographic source bins results in 55 shear-shear power spectra. Or they can be read off when running the covariance code, i.e. type 'compute_covariance_fourier 100000' and look for the output mentioning number of ggl bins accepted and/or number of cluster weka lensing bins accepted. The default numbers below most likely don't correspond to your binning choices.
nlens = 10 	# number of lens bins 
nkappa=1
nkappaxlens=10
nkappaxsource=10
nggl = 52 	# number of ggl power spectra
nshear = 55 # number of shear tomographic power spectra
ncl=15		# number of ell-bins


n2pt = (nshear+nggl+nlens)*ncl 
ndata = (nshear+nggl+nlens+nkappa+nkappaxsource+nkappaxlens)*ncl

for k in range(0,1):
  	datafile= np.genfromtxt(data[k])
  	mask = np.zeros(ndata)
	for i in range(0,datafile.shape[0]):
		if (datafile[i,1] >1.0e-15): 
			mask[i]=1.0

  	
  	covfile = np.genfromtxt(infile[k])
	cov = np.zeros((ndata,ndata))

	print ndata,n2pt,int(np.max(covfile[:,0])+1)

	for i in range(0,covfile.shape[0]):
	  	cov[int(covfile[i,0]),int(covfile[i,1])] = covfile[i,8]+covfile[i,9]
	  	cov[int(covfile[i,1]),int(covfile[i,0])] = covfile[i,8]+covfile[i,9]
		# cov[int(covfile[i,0]),int(covfile[i,1])] = covfile[i,8]
	 	# cov[int(covfile[i,1]),int(covfile[i,0])] = covfile[i,8]
	 
	
	cor = np.zeros((ndata,ndata))
	for i in range(0,ndata):
	    for j in range(0,ndata):
	    	if (cov[i,i]*cov[j,j] >0):
	       		cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])


	a = np.sort(LA.eigvals(cor[:,:]))
	print "min+max eigenvalues full cor:"
	print np.min(a), np.max(a)
	print "neg eigenvalues full cor:"
	for i in range(0,a.shape[0]):
		if (a[i]< 0.0): print a[i]


	# ############### invert shear covariance #################
	inv = LA.inv(cov[0:nshear*ncl,0:nshear*ncl])
	a = np.sort(LA.eigvals(cov[0:nshear*ncl,0:nshear*ncl]))
	print "min+max eigenvalues shear cov:"
	print np.min(a), np.max(a)
	outfile = "cov/"+outname[k]+"_shear_shear_inv"
	f = open(outfile, "w")
	for i in range(0,nshear*ncl):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,nshear*ncl):
	  		f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()

	# ############### invert kapp CMB covariance #################
	inv = LA.inv(cov[0:nkappa*ncl,0:nkappa*ncl])
	a = np.sort(LA.eigvals(cov[0:nkappa*ncl,0:nkappa*ncl]))
	print "min+max eigenvalues kappa cov:"
	print np.min(a), np.max(a)
	outfile = "cov/"+outname[k]+"_kappa_kappa_inv"
	f = open(outfile, "w")
	for i in range(0,nkappa*ncl):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,nkappa*ncl):
	  		f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()
	
	# ############### invert clustering covariance #################
	inv = LA.inv(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl])
	a = np.sort(LA.eigvals(cov[(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl,(nshear+nggl)*ncl:(nshear+nggl+nlens)*ncl]))
	print "min+max eigenvalues clustering cov:"
	print np.min(a), np.max(a)
	outfile = "cov/"+outname[k]+"_pos_pos_inv"
	f = open(outfile, "w")
	for i in range(0,nlens*ncl):
		inv[i,i]=inv[i,i]*mask[(nshear+nggl)*ncl+i]
		for j in range(0,nlens*ncl):
	  		f.write("%d %d %e\n" %(i,j, inv[i,j]))
	f.close()

	# ############### invert 2pt covariance #################
	a = np.sort(LA.eigvals(cov[0:n2pt,0:n2pt]))
	print "min+max eigenvalues 3x2pt cov:"
	print np.min(a), np.max(a)
	inv = LA.inv(cov[0:n2pt,0:n2pt])
	outfile = "cov/"+outname[k]+"_3x2pt_inv" 
	f = open(outfile, "w")
	for i in range(0,n2pt):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,n2pt):
	  		f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()

 	############### invert 6x2 covariance #################
	a = np.sort(LA.eigvals(cov[0:ndata,0:ndata]))
	print "min+max eigenvalues 6x2pt cov:"
	print np.min(a), np.max(a)
	inv = LA.inv(cov[0:ndata,0:ndata])
	outfile = "cov/"+outname[k]+"_6x2pt_inv" 
	f = open(outfile, "w")
	for i in range(0,ndata):
		inv[i,i]=inv[i,i]*mask[i]
	  	for j in range(0,ndata):
	  		f.write("%d %d %e\n" %( i,j, inv[i,j]))
	f.close()

	

	labels = (r'$C^{\gamma \gamma}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$',r'$C^{\kappa \kappa}\left(\ell,z_{\mathrm{s}_i},z_{\mathrm{s}_j}\right)$')
	ticks = np.zeros(6)
	tickx = np.zeros(5)
	ticks[1] = nshear*ncl
	ticks[2] = (nshear+nggl)*ncl
	ticks[3] = n2pt
	ticks[4] = n2pt+nkappaxlens*ncl
	ticks[5] = n2pt+(nkappaxlens+nkappaxsource)*ncl
	ticks[6] = ndata

	cor = np.zeros((ndata,ndata))
	for i in range(0,ndata):
		for j in range(0,ndata):
			if (cov[i,i]*cov[j,j] >0):
				cor[i,j] = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j])

	fs= 10
	for i in range(0,6):
  		tickx[i] = 0.5*(ticks[i]+ticks[i+1])
  		plt.plot([ticks[i]-0.5,ticks[i]-0.5],[-.5,ndata-0.5],linestyle ='-',color = 'k')
  		plt.plot([-.5,ndata-0.5],[ticks[i]-0.5,ticks[i]-0.5],linestyle ='-',color = 'k')

	plt.subplot(1, 1, 1)
	ax = plt.gca()
	im = ax.imshow(np.log10(cov[:,:]), interpolation="nearest",vmin=-25, vmax=-10)
	plt.xticks(tickx, labels,fontsize=fs)
	plt.yticks(tickx-0.5, labels,fontsize=fs)
	plt.tick_params(axis = 'x',length = 0, pad = 15)
	plt.tick_params(axis = 'y',length = 0, pad = 5)

	plt.colorbar(im)
	plt.show()
	
	print ticks
	
	# plt.figure()
	# #plt.imshow(np.log10(cov[0:1500,2000:]), interpolation="nearest",vmin=-25, vmax=-10)
	# plt.imshow(np.log10(cov[:,:]), interpolation="nearest",vmin=-25, vmax=-10)
	# #plt.imshow(cor[n2ptcl:n2ptcl+200,300:nshear*ncl], interpolation="nearest",vmax=0.5)
	# plt.colorbar()
	# plt.show()

	




