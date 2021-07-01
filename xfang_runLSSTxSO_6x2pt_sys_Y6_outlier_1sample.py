#!/usr/bin/env python

import sys
sys.path.append('/home/u1/xfang/LSSTxSO')

from cosmolike_libs_LSSTxSO_1sample import * 
from schwimmbad import MPIPool

inv=['invcov_Y1_6x2pt_1sample','invcov_Y6_6x2pt_1sample']

data=['6x2pt_LSSTxSO_Y1_1sample_dmo_outext','6x2pt_LSSTxSO_Y6_1sample_dmo_outext']

bary=['LPC_6x2pt_LSSTxSO_Y1_1sample','LPC_6x2pt_LSSTxSO_Y6_1sample']

source_z=['src_LSSTY1','src_LSSTY6'] 

lens_z=source_z # lens=src

shear_prior=[0.013,0.003] 
delta_z_prior_shear=[0.002,0.001]
delta_z_prior_clustering=delta_z_prior_shear
sigma_z_shear=[0.05,0.05]
sigma_z_clustering=sigma_z_shear
sigma_z_prior_shear=[0.006,0.003]
sigma_z_prior_clustering=sigma_z_prior_shear

nsource_table=[11.0,23.0]  
nlens_table=[18.0,41.0]
area_table=[12300.0,16500.0]

survey_designation=["LSSTxSO_Y1","LSSTxSO_Y6"]
tomo_binning_source=["source_std","source_std"]
tomo_binning_lens=["lens=src","lens=src"]

Qprior_sigma=[28.,3.5,1.7]

model=1 
file_source_z = os.path.join(dirname, "zdistris/",source_z[model])
file_lens_z = os.path.join(dirname, "zdistris/",lens_z[model])
data_file = os.path.join(dirname, "datav/",data[model])
cov_file = os.path.join(dirname, "cov/",inv[model])
#cov_file = os.path.join("/Users/timeifler/Dropbox/cosmolike_store/LSST_emu/inv/",inv[model])
chain_file = os.path.join(dirname, "chains/LSSTxSO_6x2pt_outlier_model_%d_1sample" %model)
bary_file=os.path.join(dirname, "baryons/",bary[model])

initcosmo("halofit".encode('utf-8'))
initbins(15,20.0,3000.0,3000.0,21.0,10,10)
initpriors(shear_prior[model],sigma_z_shear[model],delta_z_prior_shear[model],sigma_z_prior_shear[model],sigma_z_clustering[model],delta_z_prior_clustering[model],sigma_z_prior_clustering[model],3.0,1.2,3.8,2.0,Qprior_sigma[0],Qprior_sigma[1],Qprior_sigma[2])
initsurvey(survey_designation[model].encode('utf-8'),nsource_table[model],nlens_table[model],area_table[model])
initgalaxies(file_source_z.encode('utf-8'),file_lens_z.encode('utf-8'),"outlier_model".encode('utf-8'),"outlier_model".encode('utf-8'),tomo_binning_source[model].encode('utf-8'),tomo_binning_lens[model].encode('utf-8'))
initia("NLA_HF".encode('utf-8'),"GAMA".encode('utf-8'))

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("6x2pt".encode('utf-8'))
initdatainv(cov_file.encode('utf-8') ,data_file.encode('utf-8'), bary_file.encode('utf-8'))
initcmb("so_Y5".encode('utf-8'))
#sample_params= sample_cosmology_only()
sample_params = sample_cosmology_3x2_allsys_1sample(get_N_tomo_shear())

Nwalker = int(sys.argv[1])
sample_main_1sample(sample_params,sigma_z_shear[model],8000,Nwalker,1,chain_file, blind=False, pool=MPIPool())

