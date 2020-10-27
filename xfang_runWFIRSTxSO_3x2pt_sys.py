#!/usr/bin/env python

import sys
sys.path.append('/home/u1/xfang/LSSTxSO')

from cosmolike_libs_WFIRSTxSO import * 
from schwimmbad import MPIPool

inv=['invcov_wfirstwidexso_3x2pt','invcov_Y6_3x2pt']

data=['3x2pt_WFIRSTwidexSO5_dmo','3x2pt_LSSTxSO_Y6_dmo']

bary=['LPC_3x2pt_WFIRSTwidexSO5','LPC_3x2pt_LSSTxSO_Y6']

source_z=['zdistri_WFIRST_LSST_lensing_fine_bin','src_LSSTY6'] 

lens_z=['zdistri_WFIRST_LSST_clustering_fine_bin','lens_LSSTY6']

shear_prior=[0.002,0.002] 
delta_z_prior_shear=[0.001,0.001]
delta_z_prior_clustering=[0.001,0.001]
sigma_z_shear=[0.02,0.01]
sigma_z_clustering=[0.02,0.01]
sigma_z_prior_shear=[0.002,0.002]
sigma_z_prior_clustering=[0.002,0.003]

nsource_table=[48.0,51.0]  
nlens_table=[52.0,66.0]
area_table=[18000.0,2000.0]

survey_designation=["WFIRSTwidexSO","WFIRSTstdxSO"]
tomo_binning_source=["source_std","source_std"]
tomo_binning_lens=["WF_SN10","WF_SN10"]

Qprior_sigma=[37.,4.7,2.8]

model=0
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model])
cov_file = os.path.join(dirname, "cov/",inv[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
chain_file = os.path.join(dirname, "chains/%s_6x2pt_model" %survey_designation[model])
bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit")
initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initpriors(shear_prior[model],sigma_z_shear[model],delta_z_prior_shear[model],sigma_z_prior_shear[model],sigma_z_clustering[model],delta_z_prior_clustering[model],sigma_z_prior_clustering[model],3.0,1.2,3.8,2.0,Qprior_sigma[0],Qprior_sigma[1],Qprior_sigma[2]);
initsurvey(survey_designation[model],nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian",tomo_binning_source[model],tomo_binning_lens[model])
initia("NLA_HF","GAMA")

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("3x2pt")
initdatainv(cov_file ,data_file, bary_file)
# initcmb("so_Y5")
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_3x2_allsys(get_N_tomo_shear(),get_N_tomo_clustering())

sample_main(sample_params,sigma_z_shear[model],sigma_z_clustering[model],8000,1120,1,chain_file, blind=False, pool=MPIPool())

