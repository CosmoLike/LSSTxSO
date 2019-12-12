#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

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
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/covariances_3D.c"
#include "../cosmolike_core/theory/covariances_fourier.c"
#include "../cosmolike_core/theory/covariances_cluster.c"
#include "init_emu.c"

void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start);
void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start);
void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start);
void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster);
void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);


void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start)
{
  int nN1, nN2,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      i += Cluster.N200_Nbin*nzc1+nN1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      cov =cov_N_N(nzc1,nN1, nzc2, nN2);
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,0.0,0.0, nzc1, nN1, nzc2, nN2,cov,0.0);
    }
  }
  fclose(F1);
}

void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
     for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
       i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
       i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
       j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
       j += Cluster.N200_Nbin*nzc2+nN2;

       cov =cov_cgl_N(ell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2);
       fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell_Cluster[nl1], 0., nzc1, nzs1, nzc2, nN2,cov,0.);
     }
   }
 }
 fclose(F1);
}

void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  nzc2 = ZC(N2);
  nzs2 = ZSC(N2);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
      for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
        for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
          i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
          j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

          c_g = 0;
          c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
          if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
          fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nzs1, nzc2, nzs2,c_g, c_ng);
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j, N1,N2;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s",PATH,OUTFILE);
  F1 =fopen(filename,"w");
  for (N1 = 0; N1 < tomo.cgl_Npowerspectra; N1 ++){
    for (N2 = 0; N2 < tomo.cgl_Npowerspectra; N2 ++){
      nzc1 = ZC(N1);
      nzs1 = ZSC(N1);
      nzc2 = ZC(N2);
      nzs2 = ZSC(N2);
      for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
        for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
          for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
            for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
              i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
              j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;         
              c_g = 0;
              c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
              if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
              fprintf(F1,"%d %d %e %e %d %d %d  %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nN1,nzs1, nzc2, nN2,nzs2,c_g, c_ng);
            }
          }
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int nz1,nz2, nN2, nl1, nzc1, i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nz1 = Z1(N1);
  nz2 = Z2(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      cov = 0.;
      i = like.Ncl*N1+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      if (ell[nl1] < like.lmax_shear){cov =cov_shear_N(ell[nl1],nz1,nz2, nzc2, nN2);}
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., nz1, nz2, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN1, nN2, nzs1, nzs2,nl2, nzc2, nzs3,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzs1 = Z1(N1);
  nzs2 = Z2(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*N1+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;          
        c_g = 0.;
        c_ng = 0.;
        if (ell[nl1] < like.lmax_shear){
          c_ng = cov_NG_shear_cgl(ell[nl1],ell_Cluster[nl2],nzs1, nzs2, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.001){ 
            c_g =cov_G_shear_cgl(ell[nl1],dell_Cluster[nl2],nzs1,nzs2, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], nzs1, nzs2, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        cov =cov_ggl_N(ell[nl1],zl,zs, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., zl, zs, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2, zl, zs, nzs1, nl2, nzc2, nzs3,i,j;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],zl);
        if (weight){
          c_ng = cov_NG_ggl_cgl(ell[nl1],ell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_ggl_cgl(ell[nl1],dell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);}
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], zl, zs, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell,int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],N1);
      if (weight){
        cov =cov_cl_N(ell[nl1],N1,N1,nzc2,nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., N1, N1, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2,nzc2, nzs3,i,j,nl2;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],N1);
        if (weight){
          c_ng = cov_NG_cl_cgl(ell[nl1],ell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){
            c_g =cov_G_cl_cgl(ell[nl1],dell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          }
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
        //printf("%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
      }
    }
  }
  fclose(F1);
}

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight;
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
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
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
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,weight;
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
      if (z1 == zl){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
        if (weight){
          c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);
          if (nl1 == nl2){
            c_g =  cov_G_cl_gl_tomo(ell[nl1],dell[nl1],z1,z2,zl,zs);
          }
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
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
      if (z1 == z3){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
        if (weight) {
          c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_cl_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2, weight;
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
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
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
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n",like.Ncl*n1+nl1,like.Ncl*(n2)+nl2,ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

/*** CMBkappa part ****/
// gs_kk
void run_cov_gs_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int start)
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
void run_cov_gs_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int zs2, int start)
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
// gk_gk
void run_cov_gk_gk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zl2, int start)
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
void run_cov_gk_gs(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2, int start)
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
void run_cov_gk_kk(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int start)
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
void run_cov_gk_ks(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int zs2, int start)
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
void run_cov_gk_ss(char *OUTFILE, char *PATH, double *ell, double *dell, int zl1, int n2,int start)
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

int main(int argc, char** argv)
{
  
  int i,l,m,n,o,s,p,nl1,t,k;
  char OUTFILE[400],filename[400],arg1[400],arg2[400];
  
  int N_scenarios=3;
  double area_table[3]={2000.0,18000.0,18000.0};
  double nsource_table[3]={51.0,27.0,43.0};
  
  double nlens_table[3]={66.0,48.0,50.0};
  
  char survey_designation[3][200]={"WFIRST","LSST_Y10","WFIRST"};
  char source_zfile[3][400]={"zdistri_WFIRST_LSST_lensing_fine_bin","SRD_zdistri_model_z0=1.100000e-01_beta=6.800000e-01_Y10_source.txt","zdistri_WFIRST_LSST_lensing_fine_bin"};

  char lens_zfile[3][400]={"zdistri_WFIRST_LSST_clustering_fine_bin","SRD_zdistri_model_z0=2.800000e-01_beta=9.000000e-01_Y10_lens.txt","zdistri_WFIRST_LSST_clustering_fine_bin"};

// WFIRST when using SNR=5 clustering sample:
// double scenario_table[1][3]={{2000.0,51.0,107.0}};
// when using SNR=10 clustering sample
// double scenario_table[1][3]={{2000.0,51.0,66.0}};
// when using SNR=10 clustering sample LSST area and f_LSST of 73% in the clustering SN>5 sample, 79% 
  
   
  int hit=atoi(argv[1]);
  Ntable.N_a=20;
  k=1;
   
  for(t=2;t<3;t++){
   
    //RUN MODE setup
    init_cosmo_runmode("emu");
    init_binning_fourier(15,20.0,3000.0,3000.0,21.0,10,10);
    init_survey(survey_designation[t],nsource_table[t],nlens_table[t],area_table[t]);
    sprintf(arg1,"zdistris/%s",source_zfile[t]);
    sprintf(arg2,"zdistris/%s",lens_zfile[t]); 
    init_galaxies(arg1,arg2,"none","none",survey_designation[t]);
    init_clusters();
    init_IA("none", "GAMA"); 
    init_probes("3x2pt");
   

    //set l-bins for shear, ggl, clustering, clusterWL
    double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
    double *ell, *dell, *ell_Cluster, *dell_Cluster;
    ell=create_double_vector(0,like.Ncl-1);
    dell=create_double_vector(0,like.Ncl-1);
    ell_Cluster=create_double_vector(0,Cluster.lbin-1);
    dell_Cluster=create_double_vector(0,Cluster.lbin-1);
    int j=0;
    for(i=0;i<like.Ncl;i++){
      ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
      dell[i]=exp(log(like.lmin)+(i+1)*logdl)-exp(log(like.lmin)+(i*logdl));
      if(ell[i]<like.lmax_shear) printf("%le\n",ell[i]);
      if(ell[i]>like.lmax_shear){
        ell_Cluster[j]=ell[i];
        dell_Cluster[j]=dell[i];
        printf("%le %le\n",ell[i],ell_Cluster[j]);
        j++;
      }
    } 


    printf("----------------------------------\n");  
    
    
    sprintf(survey.name,"%s_area%le_ng%le_nl%le",survey_designation[0],survey.area,survey.n_gal,survey.n_lens);
    printf("area: %le n_source: %le n_lens: %le\n",survey.area,survey.n_gal,survey.n_lens);

//    sprintf(covparams.outdir,"/home/u17/timeifler/covparallel/"); 
    sprintf(covparams.outdir,"/halo_nobackup/cosmos/teifler/covparallel/");

    printf("----------------------------------\n");  
    sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){ 
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_shear_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
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
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
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
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
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
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
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
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_clustering_ggl(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
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
          if (fopen(filename, "r") != NULL){exit(1);}
          else {
            run_cov_ggl_shear(OUTFILE,covparams.outdir,ell,dell,l,m,k);
          }
        }
        k=k+1;
        //printf("%d\n",k);
      }
    }

    //****************************** 
    //******cluster covariance****** 
    //******************************
    
   //  sprintf(OUTFILE,"%s_nn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cluster_Nbin; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_N_N (OUTFILE,covparams.outdir,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //      //printf("%d\n",k);
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_cscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cgl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_cgl_cgl (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    } 
   //  }
   //  sprintf(OUTFILE,"%s_csn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.cgl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //          run_cov_cgl_N (OUTFILE,covparams.outdir,ell_Cluster,dell_Cluster,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    }
   //  }
   //  //shear X cluster
   //  sprintf(OUTFILE,"%s_ssn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.shear_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_shear_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
   //       }
   //     }
   //     k=k+1;
   //   }
   // } 
   // sprintf(OUTFILE,"%s_sscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   // for (l=0;l<tomo.shear_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_shear_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }
   //        k=k+1;
   //      //}
   //    }
   //  }
   //  // ggl X cluster
   //  sprintf(OUTFILE,"%s_lsn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.ggl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_ggl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);
   //        }
   //      }
   //      k=k+1;
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_lscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.ggl_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_ggl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }       
   //        k=k+1;
   //      //}
   //    }
   //  }
   //  // clustering X cluster
   //  sprintf(OUTFILE,"%s_lln_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.clustering_Npowerspectra; l++){
   //    for (m=0;m<tomo.cluster_Nbin; m++){
   //      if(k==hit){
   //        sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //        if (fopen(filename, "r") != NULL){exit(1);}
   //        else {
   //         run_cov_cl_N (OUTFILE,covparams.outdir,ell,dell,l,m,k);          
   //        }
   //      } 
   //      //printf("%d\n",k);
   //      k=k+1;
   //    }
   //  }
   //  sprintf(OUTFILE,"%s_llcs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
   //  for (l=0;l<tomo.clustering_Npowerspectra; l++){
   //    for (m=0;m<tomo.cgl_Npowerspectra; m++){
   //      //for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
   //        if(k==hit){
   //          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
   //          if (fopen(filename, "r") != NULL){exit(1);}
   //          else {
   //            run_cov_cl_cgl (OUTFILE,covparams.outdir,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
   //          }
   //        }
   //        k=k+1;
   //      //}
   //    }
   //  }
  }
  printf("number of cov blocks for parallelization: %d\n",k-1); 
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

