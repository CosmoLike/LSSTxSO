#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

#include <fftw3.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/CMBxLSS_fourier.c"
#include "../cosmolike_core/theory/covariances_CMBxLSS_fourier.c"

#include "../cosmolike_core/theory/covariances_binned_simple.c"
#include "../cosmolike_core/theory/run_covariances_fourier_binned_6x2pt.c"

#include "init_LSSxCMB.c"


int main(int argc, char** argv)
{
  
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400],arg1[400],arg2[400];
  
  int N_scenarios=2;
  double area_table[2]={12300.0,16500.0}; // Y1 corresponds to DESC SRD Y1, Y6 corresponds to assuming that we cover the full SO area=0.4*fsky and at a depth of 26.1 which is in a range of reasonable scenarios (see https://github.com/LSSTDESC/ObsStrat/tree/static/static )
  double nsource_table[2]={11.0,23.0};
  double nlens_table[2]={18.0,41.0};
  
  char survey_designation[2][200]={"LSSTxSO_Y1","LSSTxSO_Y6"};
  
  char source_zfile[2][400]={"src_LSSTY1","src_LSSTY6"};

#ifdef ONESAMPLE
  char lens_zfile[2][400]={"src_LSSTY1","src_LSSTY6"};
  nlens_table[0] = nsource_table[0];
  nlens_table[1] = nsource_table[1];
#else
  char lens_zfile[2][400]={"lens_LSSTY1","lens_LSSTY6"};
#endif

  int hit=atoi(argv[1]);
  Ntable.N_a=100;
  k=1;
  
  t = atoi(argv[2]);
  
  //RUN MODE setup
  init_cosmo_runmode("emu");
  // init_binning_fourier(20,30.0,3000.0,3000.0,21.0,10,10);
  init_binning_fourier(15,20.0,3000.0,3000.0,0.0,10,10);
  init_survey(survey_designation[t],nsource_table[t],nlens_table[t],area_table[t]);
  sprintf(arg1,"zdistris/%s",source_zfile[t]);
  sprintf(arg2,"zdistris/%s",lens_zfile[t]); 
  init_galaxies(arg1,arg2,"none","none","source_std","LSST_gold");
  init_IA("none","GAMA"); 
  init_probes("6x2pt");

  if(t==0) init_cmb("so_Y1");
  if(t==1) init_cmb("so_Y5");
  cmb.fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  //set l-bins for shear, ggl, clustering, clusterWL

  double lmin=like.lmin;
  double lmax=like.lmax;
  int Nell=like.Ncl;
  int Ncl = Nell;
  double logdl=(log(lmax)-log(lmin))/Nell;
  double *ellmin, *dell;
  ellmin=create_double_vector(0,Nell);
  dell=create_double_vector(0,Nell-1);
  double ellmax;
  for(i=0; i<Nell ; i++){
    ellmin[i]=exp(log(lmin)+(i+0.0)*logdl);
    ellmax = exp(log(lmin)+(i+1.0)*logdl);
    dell[i]=ellmax-ellmin[i];
  }
  ellmin[Nell] = ellmax;
  like.ell = ellmin;

  covparams.ng = 1;
  covparams.cng= 1;

  printf("----------------------------------\n");  
  sprintf(survey.name,"%s_area%le_ng%le_nl%le",survey_designation[t],survey.area,survey.n_gal,survey.n_lens);
  printf("area: %le n_source: %le n_lens: %le\n",survey.area,survey.n_gal,survey.n_lens);

  // sprintf(covparams.outdir,"/home/u17/timeifler/covparallel/"); 
#ifdef ONESAMPLE
  sprintf(covparams.outdir,"out_cov_lsstxso_1sample/");
#else
  sprintf(covparams.outdir,"out_cov_lsstxso/");
  //sprintf(covparams.outdir,"/halo_nobackup/cosmos/teifler/covparallel/");
#endif

  printf("----------------------------------\n");  
  if (like.shear_shear)
  {
    sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);

    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);
        }

        k=k+1;
      }
    }
  }
  if (like.pos_pos)
  {
    sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){ 
      for (m=l;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k); 
        }        
        k=k+1;
      }
    }
  }
  if (like.shear_pos)
  {
    sprintf(OUTFILE,"%s_lsls_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=l;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ggl_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }        
        k=k+1;
      }
    }
  }
  if (like.shear_pos && like.shear_shear)
  {
    sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ggl_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        } 
        k=k+1;
      }
    }
  }
  if (like.pos_pos && like.shear_shear)
  {
    sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        } 
        k=k+1;
      }
    }
  }
  if (like.pos_pos && like.shear_pos)
  {
    sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Npowerspectra; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_clustering_ggl_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.gk && like.shear_shear)
  {
    sprintf(OUTFILE,"%s_lkss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        } 
        k=k+1;
      }
    }
  }
  if (like.ks && like.shear_shear)
  {
    sprintf(OUTFILE,"%s_ksss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        } 
        k=k+1;
      }
    }
  }
  if (like.gk && like.shear_pos)
  {
    sprintf(OUTFILE,"%s_lkls_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_ggl_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.ks && like.shear_pos)
  {
    sprintf(OUTFILE,"%s_ksls_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_ggl_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.gk && like.pos_pos)
  {
    sprintf(OUTFILE,"%s_lkll_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=0;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_clustering_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.ks && like.pos_pos)
  {
    sprintf(OUTFILE,"%s_ksll_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.clustering_Npowerspectra; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_clustering_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.gk)
  {
    sprintf(OUTFILE,"%s_lklk_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.clustering_Nbin; l++){
      for (m=l;m<tomo.clustering_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_gk_gk_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.ks && like.gk)
  {
    sprintf(OUTFILE,"%s_kslk_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=0;m<tomo.clustering_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_gk_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  if (like.ks)
  {
    sprintf(OUTFILE,"%s_ksks_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Nbin; l++){
      for (m=l;m<tomo.shear_Nbin; m++){
        if(k==hit) {
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_ks_ks_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,l,m,k);  
        }
        k=k+1;
      }
    }
  }
  ////////

  if (like.kk && like.shear_shear)
  {
    sprintf(OUTFILE,"%s_kkss_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_shear_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,m,k);  
      } 
      k=k+1;
    }
  }
  if (like.kk && like.shear_pos)
  {
    sprintf(OUTFILE,"%s_kkls_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (m=0;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_ggl_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,m,k);  
      }
      k=k+1;
    }
  }
  if (like.kk && like.pos_pos)
  {
    sprintf(OUTFILE,"%s_kkll_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (m=0;m<tomo.clustering_Npowerspectra; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_clustering_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,m,k);  
      }
      k=k+1;
    }
  }
  if (like.kk && like.gk)
  {
    sprintf(OUTFILE,"%s_kklk_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (m=0;m<tomo.clustering_Nbin; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_gk_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,m,k);  
      }
      k=k+1;
    }
  }
  if (like.kk && like.ks)
  {
    sprintf(OUTFILE,"%s_kkks_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    for (m=0;m<tomo.shear_Nbin; m++){
      if(k==hit) {
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        // if (fopen(filename, "r") != NULL){exit(1);}
        run_cov_kk_ks_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,m,k);  
      }
      k=k+1;
    }
  }
  if (like.kk)
  {
    sprintf(OUTFILE,"%s_kkkk_cov_Ncl%d_Ntomo%d",survey.name,Ncl,tomo.shear_Nbin);
    if(k==hit) {
      sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
      // if (fopen(filename, "r") != NULL){exit(1);}
      run_cov_kk_kk_fourier_bin(OUTFILE,covparams.outdir,ellmin,Ncl,k);  
    }
    k=k+1;
  }


  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

