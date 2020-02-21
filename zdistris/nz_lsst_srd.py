import numpy as np 
import matplotlib.pyplot as plt 
from scipy.special import gamma

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

i=1
year = years[i]
print 'Year', year
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

Ntomo = 5

# redshift cut
nzLSS[np.logical_or(zmid<0.2, zmid>1.2)] = 0.
nzWL[zmid<0.2] = 0.

# norm_WL = WLz0**3 * gamma(3./WLalpha) / WLalpha
# norm_LSS= LSSz0**3 * gamma(3./LSSalpha) / LSSalpha

print('LSS n_eff: %f; WL n_eff: %f'%(LSSneff, WLneff))
print('In each z-bin LSS n_eff: %f; WL n_eff: %f'%(LSSneff/Ntomo, WLneff/Ntomo))

# plt.plot(zmid,nzWL, label='WL')
# plt.plot(zmid,nzLSS, label='LSS')
# plt.legend()
# plt.show()

cdf_WL = np.cumsum(nzWL)
cdf_WL /= cdf_WL[-1]
cdf_LSS = np.cumsum(nzLSS)
cdf_LSS /= cdf_LSS[-1]

nzWL1 = np.zeros(zbin)
nzWL2 = np.zeros(zbin)
nzWL3 = np.zeros(zbin)
nzWL4 = np.zeros(zbin)
nzWL5 = np.zeros(zbin)
nzLSS1 = np.zeros(zbin)
nzLSS2 = np.zeros(zbin)
nzLSS3 = np.zeros(zbin)
nzLSS4 = np.zeros(zbin)
nzLSS5 = np.zeros(zbin)

nzWL1[cdf_WL<1./Ntomo] = nzWL[cdf_WL<1./Ntomo]
nzWL2[np.logical_and(cdf_WL>=1./Ntomo, cdf_WL<2./Ntomo)] = nzWL[np.logical_and(cdf_WL>=1./Ntomo, cdf_WL<2./Ntomo)]
# print(nzWL2)
nzWL3[np.logical_and(2./Ntomo<=cdf_WL, cdf_WL<3./Ntomo)] = nzWL[np.logical_and(2./Ntomo<=cdf_WL, cdf_WL<3./Ntomo)]
nzWL4[np.logical_and(3./Ntomo<=cdf_WL, cdf_WL<4./Ntomo)] = nzWL[np.logical_and(3./Ntomo<=cdf_WL, cdf_WL<4./Ntomo)]
nzWL5[cdf_WL>=4./Ntomo] = nzWL[cdf_WL>=4./Ntomo]

nzLSS1[cdf_LSS<1./Ntomo] = nzLSS[cdf_LSS<1./Ntomo]
nzLSS2[np.logical_and(1./Ntomo<=cdf_LSS, cdf_LSS<2./Ntomo)] = nzLSS[np.logical_and(1./Ntomo<=cdf_LSS, cdf_LSS<2./Ntomo)]
nzLSS3[np.logical_and(2./Ntomo<=cdf_LSS, cdf_LSS<3./Ntomo)] = nzLSS[np.logical_and(2./Ntomo<=cdf_LSS, cdf_LSS<3./Ntomo)]
nzLSS4[np.logical_and(3./Ntomo<=cdf_LSS, cdf_LSS<4./Ntomo)] = nzLSS[np.logical_and(3./Ntomo<=cdf_LSS, cdf_LSS<4./Ntomo)]
nzLSS5[cdf_LSS>=4./Ntomo] = nzLSS[cdf_LSS>=4./Ntomo]


# plt.plot(zmid,cdf_WL, label='WL')
# plt.plot(zmid,cdf_LSS, label='LSS')
# plt.legend()
# plt.show()

Nsample = 1000000

values = np.random.rand(Nsample)
tomolist = [nzWL1,nzWL2,nzWL3,nzWL4,nzWL5,nzLSS1,nzLSS2,nzLSS3,nzLSS4,nzLSS5]

pz_tomo = np.zeros((zbin, 2*Ntomo))
for i in range(2*Ntomo):
	tomo = tomolist[i]
	cdf_tomo = np.cumsum(tomo)
	cdf_tomo /= cdf_tomo[-1]
	value_bins_tomo = np.searchsorted(cdf_tomo, values)
	z_random_from_cdf_tomo = zmid[value_bins_tomo]
	if(i<5):
		sigmaz_tomo = (1.+z_random_from_cdf_tomo)*0.05
	else:
		sigmaz_tomo = (1.+z_random_from_cdf_tomo)*0.03

	delta_z_tomo = np.random.normal(0,1,Nsample)*sigmaz_tomo

	z_convolved_tomo = z_random_from_cdf_tomo + delta_z_tomo

	## get histogram
	hist_tomo, _ = np.histogram(z_convolved_tomo, bins=zrange)


	## get normalized PDF
	pz_tomo[:,i] = hist_tomo/dz/np.sum(hist_tomo)

	z_ave = np.sum(zmid * pz_tomo[:,i]*dz)
	if(i<5):
		print('src z_ave[%d] = %f'%(i,z_ave))
		if(i==0):
			plt.plot(zmid,pz_tomo[:,i], color='C0', label='source')
		else:
			plt.plot(zmid,pz_tomo[:,i], color='C0')
	else:
		print('lens z_ave[%d] = %f'%(i-5,z_ave))
		if(i==5):
			plt.plot(zmid,pz_tomo[:,i], color='C1', label='lens')
		else:
			plt.plot(zmid,pz_tomo[:,i], color='C1')
plt.xlabel('z')
plt.ylabel('Unnormalized n(z)')
plt.legend()
plt.savefig('nz_srdLSSTY'+str(year)+'.pdf')
plt.show()

np.savetxt('source_srdLSSTY'+str(year)+'.nz', np.c_[zmid, pz_tomo[:,:5]])
np.savetxt('lens_srdLSSTY'+str(year)+'.nz', np.c_[zmid, pz_tomo[:,5:]])
