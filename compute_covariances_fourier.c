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

#include "init_LSSxCMB.c"

// Naming convention:
// l = galaxy positions ("l" as in "lens sample")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "source sample")
// And alphabetical order


void run_cov_ls_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ll_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ll_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ll_ll(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ls_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ss_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);

void run_cov_ls_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int start);
void run_cov_ls_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zs2, int start);
void run_cov_lk_lk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zl2, int start);
void run_cov_lk_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2, int start);
void run_cov_lk_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int start);
void run_cov_lk_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zs2, int start);
void run_cov_lk_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2,int start);
void run_cov_kk_kk(char *OUTFILE, char *PATH, double *ell, double *dell,int start);
void run_cov_kk_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zs2, int start);
void run_cov_kk_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n2, int start);
void run_cov_ks_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zs1, int zs2,int start);
void run_cov_ks_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int zs1, int n2, int start);
void run_cov_ll_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int start);
void run_cov_ll_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zs2, int start);
void run_cov_ll_lk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zl2,int start);


void run_cov_ls_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight,i,j;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(zl,z3)*test_zoverlap(zl,z4)){
          c_ng = cov_NG_gl_shear_tomo(ell[nl1],ell[nl2],zl,zs,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_gl_shear_tomo(ell[nl1],dell[nl1],zl,zs,z3,z4);
        }
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1;
      j=like.Ncl*n2+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ll_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4)){
          c_ng = cov_NG_cl_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
      j=like.Ncl*n2+nl2;
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ll_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
      if (weight){
        if(z1 == zl){c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);}
        if (nl1 == nl2){
          c_g =  cov_G_cl_gl_tomo(ell[nl1],dell[nl1],z1,z2,zl,zs);
        }
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
      j=like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ll_ll(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_2 = %d\n", n2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
      if ((z1 == z3) && weight) {
        c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_cl_cl_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
      j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n2)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ls_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
      if (weight && zl1 == zl2) {
        c_ng = cov_NG_gl_gl_tomo(ell[nl1],ell[nl2],zl1,zs1,zl2,zs2);
      }
      if (nl1 == nl2){
        c_g =  cov_G_gl_gl_tomo(ell[nl1],dell[nl1],zl1,zs1,zl2,zs2);
      }
      if (weight ==0 && n2 != n1){
        c_g = 0;
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1;
      j=like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2;
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_ss_ss(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear = %d\n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){
        c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_shear_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;} 
      } 
      i=like.Ncl*n1+nl1;
      j=like.Ncl*n2+nl2;     
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

/*** CMBkappa part ****/
// ls_kk
void run_cov_ls_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int start)
{
   int zl1, zs1,i,j, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for ls: %d (%d, %d)\n", n1, zl1, zs1);
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1], zl1);
         if (weight && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_gs_kk(ell[nl1],ell[nl2], zl1, zs1);
            if (nl1==nl2){
               c_g = cov_G_gs_kk(ell[nl1], dell[nl1], zl1, zs1);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,zs1,0,0,c_g,c_ng);
      }
   }
   fclose(F1);
}

// ls_ks
void run_cov_ls_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zs2, int start)
{
   int zl1, zs1,i,j,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for ls: %d (%d, %d)\n", n1, zl1, zs1);
   printf("Bin for ks: %d\n", zs2);
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1], zl1);
         if (weight && ell[nl2]<like.lmax_kappacmb) {
            if (test_zoverlap(zl1, zs2)) {
               c_ng = cov_NG_gs_ks(ell[nl1],ell[nl2], zl1, zs1, zs2);
            }
            if (nl1==nl2){
               c_g = cov_G_gs_ks(ell[nl1], dell[nl1], zl1, zs1, zs2);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,zs1,0,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}
// lk_lk
void run_cov_lk_lk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zl2, int start)
{
   int weight,i,j;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Lens bins: %d, %d\n", zl1, zl2);
   
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
         bool compute = weight * (ell[nl1]<like.lmax_kappacmb) * (ell[nl2]<like.lmax_kappacmb);
         if (compute && (zl1==zl2)) {
            c_ng = cov_NG_gk_gk(ell[nl1],ell[nl2], zl1, zl2);
         }
         if (nl1==nl2){
            c_g = cov_G_gk_gk(ell[nl1], dell[nl1], zl1, zl2);
         }
         if (!compute && (zl2!=zl1)) {
            c_g = 0;
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,0,zl2,0,c_g,c_ng);
      }
   }
   fclose(F1);
}

// lk_ls
void run_cov_lk_ls(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2, int start)
{
   int zl2,zs2,nl1,nl2,i,j, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_gkgs_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Lens bin for lk: %d\n", zl1);
   zl2 = ZL(n2); zs2 = ZS(n2);
   printf("N_tomo for ls = %d (%d, %d)\n", n2,zl2,zs2);
   
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
         if (weight && ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            if (zl1==zl2) {
               c_ng = cov_NG_gk_gs(ell[nl1],ell[nl2],zl1,zl2,zs2);
            }
            if (nl1==nl2){
               c_g = cov_G_gk_gs(ell[nl1], dell[nl1], zl1,zl2,zs2);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,0,zl2,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}

// lk_kk
void run_cov_lk_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int start)
{
   int weight,i,j;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_gkkk_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Bin for lk: %d\n", zl1);
   
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1],zl1);
         if (weight && ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_gk_kk(ell[nl1],ell[nl2],zl1);
            if (nl1==nl2){
               c_g = cov_G_gk_kk(ell[nl1],dell[nl1],zl1);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,0,0,0,c_g,c_ng);
      }
   }
   fclose(F1);
}

// lk_ks
void run_cov_lk_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zs2, int start)
{
   int weight,i,j;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_gkks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Bin for lk: %d\n", zl1);
   printf("Bin for ks: %d\n", zs2);
   
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1], zl1);
         if (weight && ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            if (test_zoverlap(zl1, zs2)) {
               c_ng = cov_NG_gk_ks(ell[nl1],ell[nl2], zl1, zs2);
            }
            if (nl1==nl2){
               c_g =  cov_G_gk_ks(ell[nl1], dell[nl1], zl1, zs2);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,0,0,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}

// lk_ss
void run_cov_lk_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2,int start)
{
   int zs2,zs3,nl1,nl2,weight,i,j;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_gkss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Bin for lk: %d\n", zl1);
   zs2 = Z1(n2); zs3 = Z2(n2);
   printf("Bin for ss: %d (%d, %d)\n", n2, zs2, zs3);
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1],zl1);
         if (weight && ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_shear){
            if (test_zoverlap(zl1, zs2)*test_zoverlap(zl1, zs3)) {
               c_ng = cov_NG_gk_ss(ell[nl1],ell[nl2],zl1,zs2,zs3);
            }
            if (nl1==nl2){
               c_g = cov_G_gk_ss(ell[nl1],dell[nl1],zl1,zs2, zs3);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl1)+nl1;
         j=like.Ncl*n2+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],zl1,0,zs2,zs3,c_g,c_ng);
      }
   }
   fclose(F1);
}

// kk_kk
void run_cov_kk_kk(char *OUTFILE, char *PATH, double *ell, double *dell,int start)
{
   int nl1,nl2,i,j,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_kkkk_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb){
            c_ng = cov_NG_kk_kk(ell[nl1],ell[nl2]);
         }
         if (nl1==nl2){
            c_g = cov_G_kk_kk(ell[nl1], dell[nl1]);
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
         fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],0,0,0,0,c_g,c_ng);
//         printf("l1, l2, Cg, Cng = %le, %le, %le, %le\n", ell[nl1], ell[nl2], c_g, c_ng);
      }
   }
   fclose(F1);
}

// kk_ks
void run_cov_kk_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zs2, int start)
{
   int nl1,nl2,i,j, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_kkks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Source bin for ks: %d\n", zs2);
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_kk_ks(ell[nl1],ell[nl2],zs2);
            if (nl1==nl2){
               c_g = cov_G_kk_ks(ell[nl1],dell[nl1], zs2);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],0,0,0,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}

// kk_ss
void run_cov_kk_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int n2, int start)
{
   int z2,z3,nl1,nl2,i,j,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   z2 = Z1(n2); z3 = Z2(n2);
   printf("Bin for ss: %d (%d, %d)\n", n2, z2, z3);
   // sprintf(filename,"%scov/%s_%s_cov_kkss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_shear){
            c_ng = cov_NG_kk_ss(ell[nl1],ell[nl2],z2,z3);
            if (nl1==nl2){
               c_g = cov_G_kk_ss(ell[nl1],dell[nl1],z2,z3);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
         j=like.Ncl*n2+nl2;
         fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],0,0,z2,z3,c_g,c_ng);
      }
   }
   fclose(F1);
}

// ks_ks
void run_cov_ks_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zs1, int zs2,int start)
{
   int nl1,nl2,i,j, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   // sprintf(filename,"%scov/%s_%s_cov_ksks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   // printf("Saving to: %s\n",filename);
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   printf("Source bins for ks: %d, %d\n", zs1, zs2);
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_ks_ks(ell[nl1],ell[nl2],zs1,zs2);
         }
         if (nl1==nl2){
            c_g = cov_G_ks_ks(ell[nl1],dell[nl1], zs1, zs2);
         }
         if ((ell[nl1]>like.lmax_kappacmb || ell[nl2]>like.lmax_kappacmb) && zs2!=zs1){
            c_g = 0;
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],0,zs1,0,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}

// ks_ss
void run_cov_ks_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int zs1, int n2, int start)
{
   int z2,z3,nl1,nl2,i,j,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   z2 = Z1(n2); z3 = Z2(n2);
   printf("Source bin for ks: %d\n", zs1);
   printf("Bin for ss: %d (%d, %d)\n", n2, z2, z3);
 
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_shear){
            c_ng = cov_NG_ks_ss(ell[nl1],ell[nl2],zs1,z2,z3);
            if (nl1==nl2){
               c_g = cov_G_ks_ss(ell[nl1],dell[nl1],zs1,z2,z3);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs1)+nl1;
         j=like.Ncl*n2+nl2;
         fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],0,zs1,z2,z3,c_g,c_ng);
      }
   }
   fclose(F1);
}


// ll_kk
void run_cov_ll_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int start)
{
   int i,j, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   // zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for ll: %d \n", n1);
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1], n1);
         if (weight && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_gg_kk(ell[nl1],ell[nl2], n1, n1);
            if (nl1==nl2){
               c_g = cov_G_gg_kk(ell[nl1], dell[nl1], n1, n1);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],n1,n1,0,0,c_g,c_ng);
      }
   }
   fclose(F1);
}

// ll_ks
void run_cov_ll_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zs2, int start)
{
   int i,j,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
   F1 =fopen(filename,"w");
   // zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for ll: %d\n", n1);
   printf("Bin for ks: %d\n", zs2);
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1], n1);
         if (weight && ell[nl2]<like.lmax_kappacmb) {
            if (test_zoverlap(n1, zs2)) {
               c_ng = cov_NG_gg_ks(ell[nl1],ell[nl2], n1, n1, zs2);
            }
            if (nl1==nl2){
               c_g = cov_G_gg_ks(ell[nl1], dell[nl1], n1, n1, zs2);
            }
         }
         i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
         j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+zs2)+nl2;
         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],n1,n1,0,zs2,c_g,c_ng);
      }
   }
   fclose(F1);
}

// ll_lk
void run_cov_ll_lk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zl2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,i,j,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  // zl = ZL(n2); zs = ZS(n2);
  zl = zl2;
  printf("N_tomo_2 = %d\n", zl2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
      if (weight){
        if(z1 == zl) {c_ng = cov_NG_gg_gk(ell[nl1],ell[nl2],z1,z2,zl);}
        if (nl1 == nl2){
          c_g =  cov_G_gg_gk(ell[nl1],dell[nl1],z1,z2,zl);
        }
      }
      i=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1;
      j=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+zl)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell[nl2],z1,z2,zl,0,c_g,c_ng);
    }
  }
  fclose(F1);
}

int main(int argc, char** argv)
{
  
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400],arg1[400],arg2[400];
  
  int N_scenarios=3;
  double area_table[3]={12300.0,16500.0,18000.}; // Y1 corresponds to DESC SRD Y1, Y6 corresponds to assuming that we cover the full SO area=0.4*fsky and at a depth of 26.1 which is in a range of reasonable scenarios (see https://github.com/LSSTDESC/ObsStrat/tree/static/static )
  double nsource_table[3]={11.0,23.0,28.0};
  double nlens_table[3]={18.0,41.0,48.0};
  
  char survey_designation[3][200]={"LSSTxSO_Y1","LSSTxSO_Y6","LSSTxSO_Y10"};
  
  char source_zfile[3][400]={"src_LSSTY1","src_LSSTY6","src_LSSTY10"};

#ifdef ONESAMPLE
  char lens_zfile[3][400]={"src_LSSTY1","src_LSSTY6","src_LSSTY10"};
  nlens_table[0] = nsource_table[0];
  nlens_table[1] = nsource_table[1];
  nlens_table[2] = nsource_table[2];
#else
  char lens_zfile[3][400]={"lens_LSSTY1","lens_LSSTY6","lens_LSSTY10"};
#endif

  int hit=atoi(argv[1]);
  Ntable.N_a=100;
  k=1;
  
  t = atoi(argv[2]);
  
  //RUN MODE setup
  init_cosmo_runmode("halofit");
  // init_binning_fourier(20,30.0,3000.0,3000.0,21.0,10,10);
  init_binning_fourier(15,20.0,3000.0,3000.0,0.0,10,10);
  init_survey(survey_designation[t],nsource_table[t],nlens_table[t],area_table[t]);
  sprintf(arg1,"zdistris/%s",source_zfile[t]);
  sprintf(arg2,"zdistris/%s",lens_zfile[t]); 
#ifdef ONESAMPLE
  init_galaxies(arg1,arg2,"none","none","source_std","lens=src");
#else
  init_galaxies(arg1,arg2,"none","none","source_std","LSST_gold");
#endif
  init_IA("none","GAMA"); 
  init_probes("6x2pt");

  if(t==0) init_cmb("so_Y1");
  if(t==1) init_cmb("so_Y5");
  if(t==2) init_cmb("so_Y5");

  //set l-bins for shear, ggl, clustering, clusterWL
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  double *ell, *dell, *ell_Cluster, *dell_Cluster;
  ell=create_double_vector(0,like.Ncl-1);
  dell=create_double_vector(0,like.Ncl-1);
  int j=0;
  for(i=0;i<like.Ncl;i++){
    ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(like.lmin)+(i+1)*logdl)-exp(log(like.lmin)+(i*logdl));
    if(ell[i]<like.lmax_shear) printf("%le\n",ell[i]);
  } 

  covparams.ng = 1;
  covparams.cng = 1;

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
  sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
    for (m=l;m<tomo.shear_Npowerspectra; m++){
      if(k==hit){ 
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ss_ss(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
      k=k+1;
      //printf("%d\n",k);
    }
  }

  sprintf(OUTFILE,"%s_lsls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=l;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ls_ls(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
     //printf("%d\n",k);
      k=k+1;
    }
  }
  sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){ //auto bins only for now!
    for (m=l;m<tomo.clustering_Npowerspectra; m++){
      if(k==hit){ 
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ll_ll(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
      k=k+1;
      //printf("%d %d %d\n",l,m,k);
    }
  }
  sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ll_ss(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.ggl_Npowerspectra; m++){
      if(k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ll_ls(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      if(k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ls_ss(OUTFILE,covparams.outdir,ell,dell,l,m,k);
      }
      k=k+1;
      //printf("%d\n",k);
    }
  }
  printf("3x2pt %d\n",k);
 // lk_lk
  sprintf(OUTFILE,"%s_lklk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Nbin; l++){
     for (m=l;m<tomo.clustering_Nbin; m++){
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_lk_lk(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        k=k+1;
     }
  }
 
  printf("lklk %d\n",k);

  // lk_ls
  sprintf(OUTFILE,"%s_lkls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Nbin; l++){
     for (m=0;m<tomo.ggl_Npowerspectra; m++){
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_lk_ls(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        k=k+1;
     }
  }
  printf("lkls %d\n",k);

  // lk_kk
  sprintf(OUTFILE,"%s_lkkk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Nbin; l++){
     if (k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_lk_kk(OUTFILE,covparams.outdir,ell,dell,l,k);
     }
     // printf("%d\n",k);
     k=k+1;
  }
  printf("lkkk %d\n",k);

  // lk_ks
  sprintf(OUTFILE,"%s_lkks_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Nbin; l++){
     for (m=0;m<tomo.shear_Nbin;m++) {
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_lk_ks(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        // printf("%d\n",k);
        k=k+1;
     }
  }
  printf("lkks %d\n",k);

  // lk_ss
  sprintf(OUTFILE,"%s_lkss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Nbin; l++){
     for (m=0;m<tomo.shear_Npowerspectra; m++){
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_lk_ss(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        k=k+1;
        // printf("%d\n",k);
     }
  }
  printf("lkss %d\n",k);

  // ls_kk
  sprintf(OUTFILE,"%s_lskk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
     if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ls_kk(OUTFILE,covparams.outdir,ell,dell,l,k);
        }
    // printf("%d\n",k);
     k=k+1;
  }
  printf("lskk%d\n",k);
  
  // ls_ks
  sprintf(OUTFILE,"%s_lsks_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
     for (m=0;m<tomo.shear_Nbin;m++) {
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ls_ks(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        // printf("%d\n",k);
        k=k+1;
     }
  }
  printf("lsks %d\n",k);
  
  // kk_kk
  sprintf(OUTFILE,"%s_kkkk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  if (k==hit){
    sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
    run_cov_kk_kk(OUTFILE,covparams.outdir,ell,dell,k);
  }
  // printf("%d\n",k);
  k=k+1;
  printf("kkkk %d\n",k);

  // kk_ks
  sprintf(OUTFILE,"%s_kkks_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Nbin; l++){
     if (k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_kk_ks(OUTFILE,covparams.outdir,ell,dell,l,k);
     }
     // printf("%d\n",k);
     k=k+1;
  }
  printf("kkks %d\n",k);

  // kk_ss
  sprintf(OUTFILE,"%s_kkss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
     if (k==hit){
      sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
      run_cov_kk_ss(OUTFILE,covparams.outdir,ell,dell,l,k);
     }
     k=k+1;
     // printf("%d\n",k);
  }
  printf("kkss %d\n",k);
  
  // ks_ks
  sprintf(OUTFILE,"%s_ksks_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Nbin; l++){
     for (m=l;m<tomo.shear_Nbin; m++){
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ks_ks(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        // printf("%d\n",k);
        k=k+1;
     }
  }
  printf("ksks %d\n",k);
  
 // ks_ss
 sprintf(OUTFILE,"%s_ksss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
 for (l=0;l<tomo.shear_Nbin; l++){
   for (m=0;m<tomo.shear_Npowerspectra; m++){
     if (k==hit){
        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
        run_cov_ks_ss(OUTFILE,covparams.outdir,ell,dell,l,m,k);
     }
     k=k+1;
     // printf("%d\n",k);
    }
  }
  printf("ksss %d\n",k);

  // ll_kk
  sprintf(OUTFILE,"%s_llkk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
     if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ll_kk(OUTFILE,covparams.outdir,ell,dell,l,k);
        }
    // printf("%d\n",k);
     k=k+1;
  }
  printf("llkk %d\n",k);
  
  // ll_ks
  sprintf(OUTFILE,"%s_llks_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
     for (m=0;m<tomo.shear_Nbin;m++) {
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ll_ks(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        // printf("%d\n",k);
        k=k+1;
     }
  }
  printf("llks %d\n",k);

  // ll_lk
  sprintf(OUTFILE,"%s_lllk_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
     for (m=0;m<tomo.clustering_Nbin;m++) {
        if (k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          run_cov_ll_lk(OUTFILE,covparams.outdir,ell,dell,l,m,k);
        }
        // printf("%d\n",k);
        k=k+1;
     }
  }
  printf("lllk %d\n",k);

  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

