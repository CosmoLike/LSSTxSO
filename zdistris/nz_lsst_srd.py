import numpy as np 
import matplotlib.pyplot as plt 

zbin=300
zrange = np.linspace(0.000001,3.5,zbin+1)
dz = zrange[1]-zrange[0]
zmin=zrange[0:300]
zmax=zrange[1:301]
zmid=np.zeros(zbin)
for i in range(zbin): 
    zmid[i]=(zmin[i]+zmax[i])/2.
nzWL=np.zeros(zbin)
nzLSS=np.zeros(zbin)

years = [1,3,6,10]
idepths = np.array([25.1, 25.7, 26.1, 26.35])

for i in range(4):
	year = years[i]
	print('Year', year)
	idepth = idepths[i]
	ilim = idepth - 1.

	LSSneff=37.8*10**(0.359*(ilim - 25.))
	LSSz0 = 0.00627*(ilim-25)**2+0.0188*(ilim-25)+0.272
	LSSalpha = 0.0125*(ilim-25)**2-0.025*(ilim-25)+0.909
	print('LSS z0: %f, alpha: %f'%(LSSz0,LSSalpha))

	# WLneff = 4.33*(idepth-25)**2+7.03*(idepth-25)+10.49
	WLneff = 10.47*10**(0.3167*(idepth - 25.))
	WLz0 = -0.0125*(idepth-25)+0.193 
	WLalpha = -0.069*(idepth-25)+0.876
	print('WL z0: %f, alpha: %f'%(WLz0,WLalpha))
	nzWL = zmid**2* np.exp(-((zmid/WLz0)**WLalpha))
	nzLSS = zmid**2 *np.exp(-((zmid/LSSz0)**LSSalpha))

	# redshift cut
	# nzLSS[np.logical_or(zmid<0.2, zmid>1.2)] = 0.
	# nzWL[zmid<0.2] = 0.

	# norm_WL = WLz0**3 * gamma(3./WLalpha) / WLalpha
	# norm_LSS= LSSz0**3 * gamma(3./LSSalpha) / LSSalpha

	print('LSS n_eff: %f; WL n_eff: %f'%(LSSneff, WLneff))
	# print('In each z-bin LSS n_eff: %f; WL n_eff: %f'%(LSSneff/Ntomo, WLneff/Ntomo))

	# np.savetxt("lens_LSSTY%d"%(year), np.c_[zmin, zmid, zmax, nzLSS])
	# np.savetxt("src_LSSTY%d"%(year), np.c_[zmin, zmid, zmax, nzWL])
	print('LSST Y%d finished'%(year))