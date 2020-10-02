#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/LSSTxSO')

from cosmolike_libs_LSSTxSO import * 
from schwimmbad import MPIPool

inv=['cov_LSSTxSO_Y1_6x2pt_inv','cov_LSSTxSO_Y6_6x2pt_inv']

data=['6x2pt_LSSTxSO_Y1_dmo','6x2pt_LSSTxSO_Y6_dmo']

bary=['LPC_6x2pt_LSSTxSO_Y1','LPC_6x2pt_LSSTxSO_Y6']

source_z=['src_LSSTY1','src_LSSTY6'] 

lens_z=['lens_LSSTY1','lens_LSSTY6']

shear_prior=[0.01,0.003] 
delta_z_prior_shear=[0.002,0.001]
delta_z_prior_clustering=[0.002,0.001]
sigma_z_shear=[0.05,0.05]
sigma_z_clustering=[0.03,0.03]
sigma_z_prior_shear=[0.006,0.003]
sigma_z_prior_clustering=[0.006,0.003]

nsource_table=[11.0,23.0]  
nlens_table=[18.0,41.0]
area_table=[12300.0,16500.0]

survey_designation=["LSSTxSO_Y1","LSSTxSO_Y6"]
tomo_binning_source=["source_std","source_std"]
tomo_binning_lens=["LSST_gold","LSST_gold"]

model=1 
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model])
cov_file = os.path.join(dirname, "cov/",inv[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
chain_file = os.path.join(dirname, "/extra/timeifler/LSSTxSO/chains/LSSTxSO_6x2pt_model_%d" %model)
bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit")
initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initpriors(shear_prior[model],sigma_z_shear[model],delta_z_prior_shear[model],sigma_z_prior_shear[model],sigma_z_clustering[model],delta_z_prior_clustering[model],sigma_z_prior_clustering[model],3.0,1.2,3.8,2.0,16.0,5.0,0.8);
initsurvey(survey_designation[model],nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian",tomo_binning_source[model],tomo_binning_lens[model])
initia("NLA_HF","GAMA")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("6x2pt")
initdatainv(cov_file ,data_file, bary_file)
initcmb("so_Y5")
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_3x2_allsys(get_N_tomo_shear(),get_N_tomo_clustering())

sample_main(sample_params,sigma_z_shear[model],sigma_z_clustering[model],8000,1120,1,chain_file, blind=False, pool=MPIPool())

