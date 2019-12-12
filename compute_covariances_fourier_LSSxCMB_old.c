#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <stdbool.h>

#include <gsl/gsl_errno.h>
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

#include "../../theory/basics.c"
#include "../../theory/structs.c"
#include "../../theory/parameters.c"
#include "../../emu13/emu.c"
#include "../../theory/recompute.c"
#include "../../theory/cosmo3D.c"
#include "../../theory/redshift.c"
#include "../../theory/halo.c"
#include "../../theory/HOD.c"
#include "../../theory/cosmo2D_fourier.c"
#include "../../theory/IA.c"
#include "../../theory/cluster.c"
#include "../../theory/BAO.c"
#include "../../theory/external_prior.c"
#include "../../theory/covariances_3D.c"
#include "../../theory/covariances_fourier.c"
#include "../../theory/covariances_cluster.c"
#include "../../theory/init.c"
#include "../../theory/CMBxLSS.c"
#include "../../theory/covariances_CMBxLSS_fourier.c"


void save_cov_ggl(char *PATH, double *ell, double *dell, int n1, int n2,int start);
void save_cov_ggl_shear(char *PATH, double *ell, double *dell, int n1, int n2,int start);
void save_cov_shear_shear(char *PATH, double *ell, double *dell, int n1, int n2,int start);

void save_cov_gk_gk(char *PATH, double *ell, double *dell, int zl1, int zl2, int start);
void save_cov_gk_gs(char *PATH, double *ell, double *dell, int zl1, int n2, int start);
void save_cov_kk_kk(char *PATH, double *ell, double *dell,int start);
void save_cov_ks_ks(char *PATH, double *ell, double *dell, int zs1, int zs2,int start);
void save_cov_kk_ks(char *PATH, double *ell, double *dell, int zs1, int start);
void save_cov_ks_ss(char *PATH, double *ell, double *dell, int zs1, int n2, int start);
void save_cov_kk_ss(char *PATH, double *ell, double *dell, int n2, int start);
void save_cov_gk_ss(char *PATH, double *ell, double *dell, int zl1, int n2,int start);
void save_cov_gk_kk(char *PATH, double *ell, double *dell, int zl1, int start);
void save_cov_gk_ks(char *PATH, double *ell, double *dell, int zl1, int zs2, int start);
void save_cov_gs_ks(char *PATH, double *ell, double *dell, int n1, int zs2, int start);
void save_cov_gs_kk(char *PATH, double *ell, double *dell, int n1, int start);

//===============================================================================================

// gk_gk
void save_cov_gk_gk(char *PATH, double *ell, double *dell, int zl1, int zl2, int start)
{
   int weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gkgk_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
//         fprintf(F1,"%d %d %e %e %d %d %e %e\n", like.Ncl*zl1+nl1, like.Ncl*zl2+nl2, ell[nl1],ell[nl2],zl1,zl2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*zl1+nl1, like.Ncl*zl2+nl2,ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gk_gs
void save_cov_gk_gs(char *PATH, double *ell, double *dell, int zl1, int n2, int start)
{
   int zl2,zs2,nl1,nl2, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gkgs_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   printf("Lens bin for gk: %d\n", zl1);
   zl2 = ZL(n2); zs2 = ZS(n2);
   printf("N_tomo for gs = %d (%d, %d)\n", n2,zl2,zs2);
   
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
//         fprintf(F1,"%d %d %e %e %d %d %d %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+n2)+nl2, ell[nl1],ell[nl2],zl1,zl2,zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gk_kk
void save_cov_gk_kk(char *PATH, double *ell, double *dell, int zl1, int start)
{
   int weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gkkk_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   printf("Bin for gk: %d\n", zl1);
   
   for (int nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (int nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         weight = test_kmax(ell[nl1],zl1);
         if (weight && ell[nl1]<like.lmax_kappacmb && ell[nl2]<like.lmax_kappacmb) {
            c_ng = cov_NG_gk_kk(ell[nl1],ell[nl2], zl1);
            if (nl1==nl2){
               c_g = cov_G_gk_kk(ell[nl1], dell[nl1], zl1);
            }
         }
//         fprintf(F1,"%d %d %e %e %d %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl2, ell[nl1],ell[nl2],zl1,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gk_ks
void save_cov_gk_ks(char *PATH, double *ell, double *dell, int zl1, int zs2, int start)
{
   int weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gkks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   printf("Bin for gk: %d\n", zl1);
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
//         fprintf(F1,"%d %d %e %e %d %d %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],zl1,zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gk_ss
void save_cov_gk_ss(char *PATH, double *ell, double *dell, int zl1, int n2,int start)
{
   int zs2,zs3,nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gkss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   printf("Bin for gk: %d\n", zl1);
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
//         fprintf(F1, "%d %d %e %e %d %d %d %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],zl1,zs2,zs3,c_g,c_ng);
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*zl1+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gs_gs
void save_cov_ggl(char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
   int zl1,zl2,zs1,zs2,nl1,nl2, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gsgs_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
         if (nl1==nl2){
            c_g =  cov_G_gl_gl_tomo(ell[nl1], dell[nl1],zl1,zs1,zl2,zs2);
         }
         // hack to ignore the bin but keep cov invertible
         if (weight==0 && n2!=n1){
            c_g = 0;
         }
//         fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gs_kk
void save_cov_gs_kk(char *PATH, double *ell, double *dell, int n1, int start)
{
   int zl1, zs1, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gskk_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for gs: %d (%d, %d)\n", n1, zl1, zs1);
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
//         fprintf(F1,"%d %d %e %e %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl2, ell[nl1],ell[nl2],zl1,zs1,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gs_ks
void save_cov_gs_ks(char *PATH, double *ell, double *dell, int n1, int zs2, int start)
{
   int zl1, zs1, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gsks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   zl1 = ZL(n1); zs1 = ZS(n1);
   printf("Bin for gs: %d (%d, %d)\n", n1, zl1, zs1);
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
//         fprintf(F1,"%d %d %e %e %d %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// gs_ss
void save_cov_ggl_shear(char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
   int zl,zs,z3,z4,nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_gsss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
               c_g =  cov_G_gl_shear_tomo(ell[nl1], dell[nl1],zl,zs,z3,z4);
            }
         }
//         fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}


// kk_kk
void save_cov_kk_kk(char *PATH, double *ell, double *dell,int start)
{
   int nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_kkkk_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
//         printf("l1, l2, Cg, Cng = %le, %le, %le, %le\n", ell[nl1], ell[nl2], c_g, c_ng);
      }
   }
   fclose(F1);
}

// kk_ks
void save_cov_kk_ks(char *PATH, double *ell, double *dell, int zs2, int start)
{
   int nl1,nl2, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_kkks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
//         fprintf(F1,"%d %d %e %e %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// kk_ss
void save_cov_kk_ss(char *PATH, double *ell, double *dell, int n2, int start)
{
   int z2,z3,nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   z2 = Z1(n2); z3 = Z2(n2);
   printf("Bin for ss: %d (%d, %d)\n", n2, z2, z3);
   sprintf(filename,"%scov/%s_%s_cov_kkss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
//         fprintf(F1, "%d %d %e %e %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],z2,z3,c_g,c_ng);
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// ks_ks
void save_cov_ks_ks(char *PATH, double *ell, double *dell, int zs1, int zs2,int start)
{
   int nl1,nl2, weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   sprintf(filename,"%scov/%s_%s_cov_ksks_Nell%d_Ns%d_Ng%d_%d", PATH, survey.name, cmb.name, like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
//         fprintf(F1,"%d %d %e %e %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],zs1,zs2,c_g,c_ng);
         fprintf(F1,"%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// ks_ss
void save_cov_ks_ss(char *PATH, double *ell, double *dell, int zs1, int n2, int start)
{
   int z2,z3,nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   z2 = Z1(n2); z3 = Z2(n2);
   printf("Source bin for ks: %d\n", zs1);
   printf("Bin for ss: %d (%d, %d)\n", n2, z2, z3);
   sprintf(filename,"%scov/%s_%s_cov_ksss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
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
//         fprintf(F1, "%d %d %e %e %d %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],zs1,z2,z3,c_g,c_ng);
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+zs1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}

// ss_ss
void save_cov_shear_shear(char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
   int z1,z2,z3,z4,nl1,nl2,weight;
   double c_ng, c_g;
   FILE *F1;
   char filename[300];
   z1 = Z1(n1); z2 = Z2(n1);
   z3 = Z1(n2); z4 = Z2(n2);
   printf("N_shear = %d (%d, %d), %d (%d, %d)\n",n1,z1,z2, n2,z3,z4);
   sprintf(filename,"%scov/%s_%s_cov_ssss_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
   printf("Saving to: %s\n",filename);
   F1 =fopen(filename,"w");
   for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
      for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
         c_ng = 0.; c_g = 0.;
         if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){
            c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
         }
         if (nl1 == nl2){
            c_g =  cov_G_shear_shear_tomo(ell[nl1], dell[nl1],z1,z2,z3,z4);
         if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;}
         }
//         fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
         fprintf(F1, "%d %d %e %e %e %e\n", like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n1)+nl1, like.Ncl*(tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+n2)+nl2, ell[nl1],ell[nl2],c_g,c_ng);
      }
   }
   fclose(F1);
}


//===============================================================================================

 int main(int argc, char** argv)
 {
   int i,l,m,n,o,s,p,nl1;
   int hit=atoi(argv[1]);
   char OUTFILE[400],PATH[400];
  
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
   sprintf(PATH,"./output/test/");

   //set l-bins for shear, ggl, clustering
   double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
   double *ell, *dell;
   ell=create_double_vector(0,like.Ncl-1);
   dell=create_double_vector(0,like.Ncl-1);
   for(i=0;i<like.Ncl;i++){
     ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
     dell[i]=exp(log(like.lmin)+(i+1)*logdl) - exp(log(like.lmin)+(i*logdl));
   } 

   //===============================================================================================
    struct timeval t1, t2;
    gettimeofday(&t1, NULL);

   int k = 0;
   
    // gk_gk
    for (l=0;l<tomo.clustering_Nbin; l++){
       for (m=l;m<tomo.clustering_Nbin; m++){
          if (k==-9) save_cov_gk_gk(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }
   
    // gk_gs
    for (l=0;l<tomo.clustering_Nbin; l++){
       for (m=0;m<tomo.ggl_Npowerspectra; m++){
          if (k==-9) save_cov_gk_gs(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }
  
    // gk_kk
    for (l=0;l<tomo.clustering_Nbin; l++){
       if (k==-9) save_cov_gk_kk(PATH,ell,dell,l,k);
       //printf("%d\n",k);
       k=k+1;
    }
    
    // gk_ks
    for (l=0;l<tomo.clustering_Nbin; l++){
       for (m=0;m<tomo.shear_Nbin;m++) {
          if (k==-9) save_cov_gk_ks(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }

    // gk_ss
    for (l=0;l<tomo.clustering_Nbin; l++){
       for (m=0;m<tomo.shear_Npowerspectra; m++){
          if (k==-9) save_cov_gk_ss(PATH,ell,dell,l,m,k);
          k=k+1;
          //printf("%d\n",k);
       }
    }
   
   // gs_gs
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
       for (m=l;m<tomo.ggl_Npowerspectra; m++){
          if (k==-9) save_cov_ggl(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }
    
    // gs_kk
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
       if (k==-9) save_cov_gs_kk(PATH,ell,dell,l,k);
//       printf("%d\n",k);
       k=k+1;
    }
    
    // gs_ks
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
       for (m=0;m<tomo.shear_Nbin;m++) {
          if (k==-9) save_cov_gs_ks(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }
    
    // gs_ss
    for (l=0;l<tomo.ggl_Npowerspectra; l++){
       for (m=0;m<tomo.shear_Npowerspectra; m++){
          if (k==-9) save_cov_ggl_shear(PATH,ell,dell,l,m,k);
          k=k+1;
          //printf("%d\n",k);
       }
    }
    
    // kk_kk
    if (k==k) save_cov_kk_kk(PATH,ell,dell,k);
    //printf("%d\n",k);
    k=k+1;
    
    // kk_ks
    for (l=0;l<tomo.shear_Nbin; l++){
       if (k==-9) save_cov_kk_ks(PATH,ell,dell,l,k);
       //printf("%d\n",k);
       k=k+1;
    }

    // kk_ss
    for (l=0;l<tomo.shear_Npowerspectra; l++){
       if (k==-9) save_cov_kk_ss(PATH,ell,dell,l,k);
       k=k+1;
       //printf("%d\n",k);
    }
    
    // ks_ks
    for (l=0;l<tomo.shear_Nbin; l++){
       for (m=l;m<tomo.shear_Nbin; m++){
          if (k==-9) save_cov_ks_ks(PATH,ell,dell,l,m,k);
          //printf("%d\n",k);
          k=k+1;
       }
    }
    
   // ks_ss
   for (l=0;l<tomo.shear_Nbin; l++){
     for (m=0;m<tomo.shear_Npowerspectra; m++){
       if (k==-9) save_cov_ks_ss(PATH,ell,dell,l,m,k);
       k=k+1;
       //printf("%d\n",k);
     }
   }
   
    // ss_ss
    for (l=0;l<tomo.shear_Npowerspectra; l++){
       for (m=l;m<tomo.shear_Npowerspectra; m++){
          if (k==-9) save_cov_shear_shear(PATH,ell,dell,l,m,k);
          k=k+1;
          //printf("%d\n",k);
       }
    }

    
    
    gettimeofday(&t2, NULL);
    double duration = (t2.tv_sec - t1.tv_sec)/60.;  // in min
    printf("Took %.2f min\n", duration);
    
   printf("Number of cov blocks for parallelization: %d\n",k);
   printf("Nell=%d, Ng=%d, Ns=%d, Ngs=%d, Nss=%d\n", like.Ncl, tomo.clustering_Nbin, tomo.shear_Nbin, tomo.ggl_Npowerspectra, tomo.shear_Npowerspectra);
   
   printf("-----------------\n");
   printf("PROGRAM EXECUTED\n");
   printf("-----------------\n");
   return 0;   
 }

