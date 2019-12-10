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


#include "../cosmolike_core/theory/baryons.h"
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
#include "../cosmolike_core/theory/redmagic.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/covariances_cluster.c"
#include "init_lsst_real.c"
#include "../cosmolike_core/theory/init_baryon.c"
#include "../cosmolike_core/theory/priors_mpp.c"
#include "../cosmolike_core/theory/CMBxLSS.c"
#include "../cosmolike_core/theory/covariances_CMBxLSS_fourier.c"

//===============================================================================================



double C_shear_tomo_sys(double ell,int z1,int z2);
double C_cgl_tomo_sys(double ell_Cluster,int zl,int nN, int zs);
double C_gl_tomo_sys(double ell,int zl,int zs);

void set_data_shear(double *ell, double *data, int start);
void set_data_ggl(double *ell, double *data, int start);

void compute_data_vector(char *PATH,
                         double OMM, double S8, double NS, double W0,double WA, double OMB, double H0,
                         double B1, double B2, double B3, double B4,
                         double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1,
                         double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1,
                         double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10,
                         double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q);

double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0,
                      double B1, double B2, double B3, double B4,
                      double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1,
                      double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1,
                      double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10,
                      double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q);

double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in);

int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}


//===============================================================================================
// gk

double C_gk_sys(double ell, int zl)
{
   double C;
// MANUWARNING: implement IA
   if(like.IA!=1) C = C_gk_nointerp(ell,zl);
   if(like.IA==1) C = C_gk_nointerp(ell,zl);
   if(like.IA==2) C += 0.;
   return C;
}

void set_data_gk(double *ell, double *data, int start)
{
   for (int nz=0; nz<tomo.clustering_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         bool compute = test_kmax(ell[i],nz) * (ell[i]<like.lmax_kappacmb);
         if (compute){
            data[start+(like.Ncl*nz)+i] = C_gk_sys(ell[i],nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}


//===============================================================================================
// gs

double C_gl_tomo_sys(double ell,int zl,int zs)
{
  double C;
  if(like.IA!=1) C = C_gl_tomo_nointerp(ell,zl,zs);
  if(like.IA==1) C = C_ggl_IA(ell,zl,zs);
  if(like.IA==2) C += C_gI_lin_nointerp(ell,zl,zs);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
  return C;
}

void set_data_ggl(double *ell, double *data, int start)
{
   for (int nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
      int zl = ZL(nz);
      int zs = ZS(nz);
      for (int i = 0; i < like.Ncl; i++){
         if (test_kmax(ell[i],zl)){
            data[start+(like.Ncl*nz)+i] = C_gl_tomo_sys(ell[i],zl,zs);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      } 
   }
}

//===============================================================================================
// kk

void set_data_kk(double *ell, double *data, int start)
{
   for (int i=0; i<like.Ncl; i++){
      if (ell[i]<like.lmax_kappacmb){
//MANUWARNING: implement the IA
         data[start+i] = C_kk_nointerp(ell[i]);
      }
      else{
         data[start+i] = 0.;
      }
   }
}

//===============================================================================================
// ks

double C_ks_sys(double ell, int zs)
{
   double C;
   // MANUWARNING: implement IA
   if(like.IA!=1) C = C_ks_nointerp(ell,zs);
   if(like.IA==1) C = C_ks_nointerp(ell,zs);
   if(like.IA==2) C += 0.;
   if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
   return C;
}

void set_data_ks(double *ell, double *data, int start)
{
   for (int nz=0; nz<tomo.shear_Nbin; nz++){
      for (int i=0; i<like.Ncl; i++){
         if (ell[i]<like.lmax_kappacmb) {
            data[start+(like.Ncl*nz)+i] = C_ks_sys(ell[i],nz);
         }
         else{
            data[start+(like.Ncl*nz)+i] = 0.;
         }
      }
   }
}


//===============================================================================================
// ss

double C_shear_tomo_sys(double ell, int z1, int z2)
{
   double C;
   if(like.IA!=1) C = C_shear_tomo_nointerp(ell,z1,z2);
   if(like.IA==1) C = C_shear_tomo_nointerp(ell,z1,z2)+C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
   if(like.IA==2) C += C_II_lin_nointerp(ell,z1,z2)+C_GI_lin_nointerp(ell,z1,z2);
   if(like.shearcalib==1) C *= (1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
   return C;
}

void set_data_shear(double *ell, double *data, int start)
{
  int i,z1,z2,nz;
  double a;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < like.Ncl; i++){
      if (ell[i] < like.lmax_shear){ data[like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+nz)+i] = C_shear_tomo_sys(ell[i],z1,z2);}
      else {data[like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+nz)+i] = 0.;}
    }
  }
}


//===============================================================================================

int set_cosmology_params(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=S8;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  cosmology.h0=H0;

  if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
  if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  return 1;
}

void set_nuisance_shear_calib(double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10)
{
  nuisance.shear_calibration_m[0] = M1;
  nuisance.shear_calibration_m[1] = M2;
  nuisance.shear_calibration_m[2] = M3;
  nuisance.shear_calibration_m[3] = M4;
  nuisance.shear_calibration_m[4] = M5;
  nuisance.shear_calibration_m[5] = M6;
  nuisance.shear_calibration_m[6] = M7;
  nuisance.shear_calibration_m[7] = M8;
  nuisance.shear_calibration_m[8] = M9;
  nuisance.shear_calibration_m[9] = M10;
}

int set_nuisance_shear_photoz(double SP1,double SP2,double SP3,double SP4,double SP5,double SP6,double SP7,double SP8,double SP9,double SP10,double SPS1)
{
  int i;
  nuisance.bias_zphot_shear[0]=SP1;
  nuisance.bias_zphot_shear[1]=SP2;
  nuisance.bias_zphot_shear[2]=SP3;
  nuisance.bias_zphot_shear[3]=SP4;
  nuisance.bias_zphot_shear[4]=SP5;
  nuisance.bias_zphot_shear[5]=SP6;
  nuisance.bias_zphot_shear[6]=SP7;
  nuisance.bias_zphot_shear[7]=SP8;
  nuisance.bias_zphot_shear[8]=SP9;
  nuisance.bias_zphot_shear[9]=SP10;
  
  for (i=0;i<tomo.shear_Nbin; i++){ 
    nuisance.sigma_zphot_shear[i]=SPS1;
    if (nuisance.sigma_zphot_shear[i]<0.001) return 0;
  }
  return 1;
}

int set_nuisance_clustering_photoz(double CP1,double CP2,double CP3,double CP4,double CP5,double CP6,double CP7,double CP8,double CP9,double CP10,double CPS1)
{
   int i;
   nuisance.bias_zphot_clustering[0]=CP1;
   nuisance.bias_zphot_clustering[1]=CP2;
   nuisance.bias_zphot_clustering[2]=CP3;
   nuisance.bias_zphot_clustering[3]=CP4;
   nuisance.bias_zphot_clustering[4]=CP5;
   nuisance.bias_zphot_clustering[5]=CP6;
   nuisance.bias_zphot_clustering[6]=CP7;
   nuisance.bias_zphot_clustering[7]=CP8;
   nuisance.bias_zphot_clustering[8]=CP9;
   nuisance.bias_zphot_clustering[9]=CP10;
   
   for (i=0;i<tomo.clustering_Nbin; i++){
      nuisance.sigma_zphot_clustering[i]=CPS1;
      if (nuisance.sigma_zphot_clustering[i]<0.001) return 0;
   }
   return 1;
}

int set_nuisance_ia(double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q)
{
  nuisance.A_ia=A_ia;  
  nuisance.beta_ia=beta_ia;
  nuisance.eta_ia=eta_ia;
  nuisance.eta_ia_highz=eta_ia_highz;
  nuisance.LF_alpha=LF_alpha;
  nuisance.LF_P=LF_P;
  nuisance.LF_Q=LF_Q;
  nuisance.LF_red_alpha=LF_red_alpha;
  nuisance.LF_red_P=LF_red_P;
  nuisance.LF_red_Q=LF_red_Q;
  if (nuisance.A_ia < 0.0 || nuisance.A_ia > 10.0) return 0;
  if (nuisance.beta_ia < -1.0 || nuisance.beta_ia > 3.0) return 0;
  if (nuisance.eta_ia < -3.0 || nuisance.eta_ia> 3.0) return 0;
  if (nuisance.eta_ia_highz < -1.0 || nuisance.eta_ia_highz> 1.0) return 0;
  if(like.IA!=0){
   if (check_LF()) return 0;
  }
return 1;
}

int set_nuisance_gbias(double B1, double B2, double B3, double B4)
{
  int i;
  gbias.b[0] = B1;
  gbias.b[1] = B2;
  gbias.b[2] = B3;
  gbias.b[3] = B4;
  if(like.bias==1){
    for (i = 0; i < 4; i++){
      if (gbias.b[i] < 0.8 || gbias.b[i] > 2.0) return 0;
    }
  }
  return 1;
}


//===============================================================================================

void compute_data_vector(char *PATH,
                         double OMM, double S8, double NS, double W0,double WA, double OMB, double H0,
                         double B1, double B2, double B3, double B4,
                         double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1,
                         double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1,
                         double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10,
                         double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q)
{
   
   int i,j,k,m=0,l;
   static double *pred;
   static double *ell;
   static double darg;
   double chisqr,a,log_L_prior=0.0;
   
   if(ell==0){
      pred= create_double_vector(0, like.Ndata-1);
      ell= create_double_vector(0, like.Ncl-1);
      darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
      for (l=0;l<like.Ncl;l++){
         ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
      }
   }
   
   if (set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0)==0) printf("Cosmology out of code boundaries\n");
   set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
   set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1);
   set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1);
   set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q);
   set_nuisance_gbias(B1,B2,B3,B4);
   
   int start=0;
   if (like.gk) {
      printf("Computing data vector: gk\n");
      set_data_gk(ell, pred, start);
      start += like.Ncl * tomo.clustering_Nbin;
   }
   if (like.shear_pos) {
      printf("Computing data vector: gs\n");
      set_data_ggl(ell, pred, start);
      start += like.Ncl * tomo.ggl_Npowerspectra;
   }
   if (like.kk) {
      printf("Computing data vector: kk\n");
      set_data_kk(ell, pred, start);
      start += like.Ncl;
   }
   if (like.ks) {
      printf("Computing data vector: ks\n");
      set_data_ks(ell, pred, start);
      start += like.Ncl * tomo.shear_Nbin;
   }
   if(like.shear_shear==1) {
      printf("Computing data vector: ss\n");
      set_data_shear(ell, pred, start);
      start += like.Ncl * tomo.shear_Npowerspectra;
   }
   
   FILE *F;
   char filename[300];
   sprintf(filename,"%sdatav/%s_%s_%s_Nell%d_Ns%d_Ng%d",PATH,survey.name,cmb.name,like.probes, like.Ncl,tomo.shear_Nbin, tomo.clustering_Nbin);
   printf("Saving data vector to %s\n", filename);
   F=fopen(filename,"w");
   for (i=0;i<like.Ndata; i++){
      fprintf(F,"%d %le\n",i,pred[i]);
   }
   fclose(F);
}

//===============================================================================================

double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0,
                      double B1, double B2, double B3, double B4,
                      double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1,
                      double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1,
                      double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10,
                      double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q)
{
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
  }
  if (set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e8;
  }
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  if (set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e8; 
  }
  if (set_nuisance_gbias(B1,B2,B3,B4)==0){
    printf("Bias out of bounds\n");
    return -1.0e8;
  }
       
  // printf("like %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  if(like.Aubourg_Planck_BAO_SN==1) log_L_prior+=log_L_Planck_BAO_SN();
  if(like.SN==1) log_L_prior+=log_L_SN();
  if(like.BAO==1) log_L_prior+=log_L_BAO();
  if(like.IA!=0) log_L_prior+=log_L_ia();
  if(like.IA!=0) log_L_prior+=log_like_f_red();
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();
  if(like.clusterMobs==1) log_L_prior+=log_L_clusterMobs();
 
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  
   int start=0;
   if (like.gk) {
//      printf("Computing data vector: gk\n");
      set_data_gk(ell, pred, start);
      start += like.Ncl * tomo.clustering_Nbin;
   }
   if (like.shear_pos) {
//      printf("Computing data vector: gs\n");
      set_data_ggl(ell, pred, start);
      start += like.Ncl * tomo.ggl_Npowerspectra;
   }
   if (like.kk) {
//      printf("Computing data vector: kk\n");
      set_data_kk(ell, pred, start);
      start += like.Ncl;
   }
   if (like.ks) {
//      printf("Computing data vector: ks\n");
      set_data_ks(ell, pred, start);
      start += like.Ncl * tomo.shear_Nbin;
   }
   if(like.shear_shear==1) {
//      printf("Computing data vector: ss\n");
      set_data_shear(ell, pred, start);
      start += like.Ncl * tomo.shear_Npowerspectra;
   }

  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
    for (j=0; j<like.Ndata; j++){
      a=(pred[i]-data_read(1,i))*invcov_read(1,i,j)*(pred[j]-data_read(1,j));
      chisqr=chisqr+a;
    }
//    if (fabs(data_read(1,i)/pred[i]-1.0) >1.e-6){
//      printf("%d %le %le %le\n",i,data_read(1,i),pred[i],data_read(1,i)/pred[i]-1.);
//    }
  }
  if (chisqr<0.0){
    printf("errror: chisqr=%le < 0\n", chisqr);
  }
  if (chisqr<-1.0) exit(EXIT_FAILURE);
  
//printf("%le\n",chisqr);
   return -0.5*chisqr+log_L_prior;
}


//===============================================================================================

void init_like(void) {

   init_cosmo();
   init_binning(10.);
   
   //    init_survey("LSST");
   //    // sources, lenses, photo-z error on lens, photo-z error on sources, clustering sample
   //    init_galaxies("../../zdistris/zdistribution_LSST", "../../zdistris/zdistribution_const_comoving", "none", "gaussian", "redmagic");
   //    init_cmb("cmbs4");
   
   init_survey("HSC");
   // sources, lenses, photo-z error on lens, photo-z error on sources, clustering sample
   init_galaxies("../../zdistris/zdistribution_HSC", "../../zdistris/zdistribution_cmassdr12_0.4-0.7", "none", "gaussian", "cmass");
   init_cmb("advact");
   
   init_IA("none", "GAMA");
   // set output path
   char PATH[300];
   sprintf(PATH,"./output/test/");
   
   // MANUWARNING: are these the right priors for me?
   //init_priors("Planck_BAO_SN_Aubourg","DES_SN","PhotoBAO");
   init_priors("none", "none", "none");
   
   init_probes("LSSxCMB");
   
   compute_data_vector(PATH,
                       0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,         // cosmo params
                       2.,2.,0.,0.,                                         // galaxy bias
                       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,        // SP
                       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,        // CP
                       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,             // shear multiplicative bias
                       5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0);         // IA
   
   // read in data vector and inverse cov matrix
   char pathInvCov[500];
// MANUWARNING: using the Gaussian cov for now!!!
   sprintf(pathInvCov,"%scov/invcovG_%s_%s_Nell%d_Ns%d_Ng%d",PATH,survey.name,cmb.name, like.Ncl,tomo.shear_Nbin, tomo.clustering_Nbin);
   char pathDataVector[500];
   sprintf(pathDataVector,"%sdatav/%s_%s_%s_Nell%d_Ns%d_Ng%d",PATH,survey.name,cmb.name,like.probes, like.Ncl,tomo.shear_Nbin, tomo.clustering_Nbin);
   init_data_inv(pathInvCov, pathDataVector);
   
   return;

}
double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in)
{

    double like = log_multi_like(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, 

        in.bias[0], in.bias[1], in.bias[2], in.bias[3], 
        in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
        in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
        in.source_z_s, 
        in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
        in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
        in.lens_z_s, 
        in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
        in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
        in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,
        in.lf[0], in.lf[1], in.lf[2], in.lf[3], in.lf[4], in.lf[5]);
  
    return like;
}



//===============================================================================================

 int main(void)
{
   
   init_like();
   
   
   // test the likelihood
   struct timeval t1, t2;
   gettimeofday(&t1, NULL);
   
   double test = log_multi_like(0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,
                                2.,2.,0.,0.,
                                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,
                                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,
                                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

   gettimeofday(&t2, NULL);
   double duration = (t2.tv_sec - t1.tv_sec);  // in sec
   duration += (t2.tv_usec - t1.tv_usec)/1.e6;  // in sec
   printf("Likelihood evaluated in %.3f sec\n", duration);

   printf("Recovered likelihood of truth (should be 0) = %le\n", test);
   
  return 0;
}


