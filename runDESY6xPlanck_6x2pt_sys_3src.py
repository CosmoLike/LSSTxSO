#!/usr/bin/env python

import sys
sys.path.append('/home/u1/xfang/CosmoLike/LSSTxSO')

from cosmolike_libs_DESxPlanck import * 
from schwimmbad import MPIPool

cov='cov_desy6xplanck_3src'

data='6x2pt_desy6xplanck_6x2pt_3bins.txt'
mask='xi_Y6_6x2pt_3src.mask'
# bary=['LPC_6x2pt_LSSTxSO_Y1','LPC_6x2pt_LSSTxSO_Y6']

source_z='source_DESY6_3bins.nz'

lens_z='mcal_1101_lens.nz'

shear_prior=0.005
delta_z_prior_shear=0.002
delta_z_prior_clustering=0.005
# sigma_z_shear=[0.05,0.05]
# sigma_z_clustering=[0.03,0.03]
# sigma_z_prior_shear=[0.006,0.003]
# sigma_z_prior_clustering=[0.006,0.003]


file_source_z = os.path.join(dirname, "zdistris/",source_z)
file_lens_z = os.path.join(dirname, "zdistris/",lens_z)
data_file = os.path.join(dirname, "datav/",data)
cov_file = os.path.join(dirname, "cov/",cov)
mask_file = os.path.join(dirname, "datav/",mask)
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
chain_file = os.path.join(dirname, "chains/DESY6xPlanck_6x2pt_3src")
#chain_file = os.path.join(dirname, "/Users/timeifler/Dropbox/cosmolike_store/LSSTxSO/chains/LSSTxSO_6x2pt_model_%d" %model)
# bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit")
initbins(20,30.0,3000.0,3000.0,20.0)
initia(4)
init_source_sample(file_source_z,4)
init_lens_sample(file_lens_z,5)

initpriors(shear_prior,delta_z_prior_shear,delta_z_prior_clustering);
# initsurvey(survey_designation[model],nsource_table[model],nlens_table[model],area_table[model])
# initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian",tomo_binning_source[model],tomo_binning_lens[model])
# initia("NLA_HF","GAMA")

initprobes("6x2pt")
initdatacovmask(cov_file ,data_file, mask_file)
initcmb("planck")
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_3x2_allsys(get_N_tomo_shear(),get_N_tomo_clustering())

sample_main(sample_params,4000,560,1,chain_file, blind=False, pool=MPIPool())

