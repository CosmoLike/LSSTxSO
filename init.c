double invcov_read(int READ, int ci, int cj);
double invcov_mask(int READ, int ci, int cj);
double mask(int READ, int ci);
double data_read(int READ, int ci);
void init_data_real(char *COV_FILE, char *MASK_FILE, char *DATA_FILE);
void init_data_inv(char *INV_FILE, char *DATA_FILE);
void init_priors(char *cosmoPrior1, char *cosmoPrior2, char *cosmoPrior3, char *cosmoPrior4);
void init_survey(char *surveyname);
void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample);
void init_cosmo();
void init_cosmo_runmode(char *runmode);
void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int cluster_rich_Nbin);
void init_binning_real(int Nt, double min, double max);
void init_probes(char *probes);
void init_probes_real(char *probes);

void set_galaxies_LSST(double density);
void set_galaxies_DES(double density);
void set_galaxies_LSST_like(double density); //DESC LSST-like specfications; for DESC forecasts
void set_galaxies_LSST_gold(double density);
void set_galaxies_DES_Y1(void);
void set_galaxies_DES_SV(double density);
void set_galaxies_CMASS(double density); // MANUWARNING: not yet really implemented!
void set_galaxies_source();
void set_clusters_LSST(); //set parameters for LSST/WFIRST forecasts
void set_clusters_DES(); //set parameters for DES forecasts
void init_wlphotoz_stage3();
void init_wlphotoz_stage4();
void init_lens_sample(char *lensphotoz, char *galsample);
void init_source_sample(char *sourcephotoz);
void init_lens_sample_();
void init_source_sample_();
void init_baryons();
void init_clphotoz_redmagic();
void init_clphotoz_benchmark();
void init_clphotoz_cmass();
void init_clphotoz_LSST_gold();
void init_clphotoz_source();
void init_clusterMobs();
void set_equal_tomo_bins(int Ntomo);
void init_IA(char *model,char *lumfct);
void init_HOD_rm();
void init_Pdelta();

void init_cmb();
void set_cmb_actpol();
void set_cmb_advact();
void set_cmb_cmbs4();

double mask(int READ, int ci)
{
  int i,intspace;
  static double *mask =0;
  if(READ==0 || mask ==0){
    FILE *F;
    mask  = create_double_vector(0, like.Ndata-1); 
    double *maskc;
    maskc  = create_double_vector(0, like.Ndata-1); 
    F=fopen(like.MASK_FILE,"r");
      if (!F){
        printf("init.c: invcov_mask: like.MASK_FILE = %s not found!\nEXIT!\n",like.MASK_FILE);
        exit(1);
      }
      int N = 0;
      for (i=0;i<like.Ndata; i++){
        fscanf(F,"%d %le\n",&intspace,&mask[i]);
        maskc[i] = mask[i];
        N += mask[i];
        if(i==399) printf("WL %d bins within angular mask\n",N);
        if(i==699) printf("WL+GGL %d bins within angular mask\n",N);
      }
     fclose(F);
     printf("%d bins within angular mask\n",N);
     printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d, like.ks = %d, like.gk = %d\n\n",like.pos_pos,like.shear_pos,like.shear_shear,like.ks,like.gk); 
     int N3x2pt, N5x2pt;
     N3x2pt = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);    
     N5x2pt = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin+tomo.clustering_Nbin);
    //test whether Ndata assumes 3x2pt or 5x2pt format
    //if so, mask out probes excluded from the analysis
     if (N == N3x2pt || N== N5x2pt){
      if(like.shear_shear==0){
        printf("masking out shear-shear bins\n");
       for (i = 0; i< like.Ntheta*2*tomo.shear_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.pos_pos==0){
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra);
        printf("masking out clustering bins\n");
        for (i = N; i< N+like.Ntheta*tomo.clustering_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.shear_pos==0){
        N = like.Ntheta*2*tomo.shear_Npowerspectra;
        printf("masking out ggl bins\n");
        for (i = N; i <N+like.Ntheta*tomo.ggl_Npowerspectra; i++){mask[i] = 0.;}
      }
    }
    //test whether Ndata 5x2pt format
    //if so, mask out probes excluded from the analysis
    if (like.Ndata == N5x2pt){
      if(like.ks==0){
        printf("masking out shear x kappa bins\n");
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
        for (i = N; i <N+like.Ntheta*tomo.shear_Nbin; i++){mask[i] = 0.;}
      }
      if(like.gk==0){
        printf("masking out galaxies x kappa bins\n");
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin);
        for (i = N; i < N+like.Ntheta*tomo.clustering_Nbin; i++){mask[i] = 0.;}
      }
    }
    N = 0;
    for (i=0;i<like.Ndata; i++){
      //printf("mask(%d) = %.1f (was %.1f before probe cut)\n",i,mask[i],maskc[i]);
      N +=  mask[i];
    }
    printf("%d data points left after masking probes\n",N);
    if (N == 0){
      printf("init.c: mask: no data points left\nEXIT\n");
      exit(1);
    }
    printf("READ MASK FILE\n");
  }
  return mask[ci];
}

int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line [1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = &line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}

double invcov_mask(int READ, int ci, int cj)
{
  static double **inv =0;

  if(READ==0 || inv ==0){
    double cov_G,cov_NG,doublespace,m;
    int i,j,intspace;
    gsl_matrix * cov   = gsl_matrix_calloc(like.Ndata, like.Ndata);      
    gsl_matrix * inv_c   = gsl_matrix_calloc(like.Ndata, like.Ndata);
    gsl_matrix_set_zero(cov);
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
      
    FILE *F;
   int n_rows =count_rows(like.COV_FILE,' ');
   F=fopen(like.COV_FILE,"r");
   if (!F){printf("init.c: invcov_mask: like.COV_FILE = %s not found!\nEXIT!\n",like.COV_FILE);exit(1);}
   switch (n_rows){
    case 3: while (fscanf(F,"%d %d %le\n", &i, &j, &cov_G) ==3) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(1,i)*mask(1,j);}
            //  printf("%d %d (/%d) %e\n",i,j,like.Ndata,cov_G);
              gsl_matrix_set(cov,i,j,(cov_G)*m);
              gsl_matrix_set(cov,j,i,(cov_G)*m);
            }
          } break;
    case 4: while (fscanf(F,"%d %d %le %le\n", &i, &j, &cov_G, &cov_NG) ==4) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(1,i)*mask(1,j);}
              //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
              gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
              gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
            }
          } break;
    case 10: while (fscanf(F,"%d %d %le %le %d %d %d %d %le %le\n", &i, &j, &doublespace, &doublespace,&intspace,&intspace,&intspace,&intspace,&cov_G,&cov_NG) ==10) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(1,i)*mask(1,j);}
              //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
              gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
              gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
            }
          } break;
    default: printf("init.c:invcov_mask: covariance file %s has %d columns - unsupported format!\nEXIT\n",like.COV_FILE,n_rows);exit(1);
   }
   fclose(F);
   printf("READ COV_FILE\n");
   //SVD_inversion(cov,inv_c,like.Ndata);
  invert_matrix_colesky(cov);
   // printf("cov inverted\n");
   for (i=0;i<like.Ndata; i++){
    for (j=0;j<like.Ndata; j++){
      //apply mask again, to make sure numerical errors in matrix inversion don't cause problems...
      //also, set diagonal elements corresponding to datavector elements outside mask to zero, so that these elements don't contribute to chi2
      inv[i][j] =gsl_matrix_get(cov,i,j)*mask(1,i)*mask(1,j);
    }
   }
   gsl_matrix_free(cov);
   gsl_matrix_free(inv_c);
   printf("FINISHED BUILDING INV COV\n");
 }    
 return inv[ci][cj];
}

double invcov_read(int READ, int ci, int cj)
{
  int i,j,intspace;
  static double **inv =0;

  if(READ==0 || inv == 0){
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.INV_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      for (j=0;j<like.Ndata; j++){
       fscanf(F,"%d %d %le\n",&intspace,&intspace,&inv[i][j]);  
     }
   }
   fclose(F);
   printf("FINISHED READING COVARIANCE\n");
 }    
 return inv[ci][cj];
}


double data_read(int READ, int ci)
{
  int i,intspace;
  static double *data = 0;
  
  if(READ==0 || data ==0){
    data  = create_double_vector(0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.DATA_FILE,"r");
    for (i=0;i<like.Ndata; i++){  
      fscanf(F,"%d %le\n",&intspace,&data[i]);
    }
    fclose(F);
    printf("FINISHED READING DATA VECTOR\n");
  }    
  return data[ci];
}


void init_cosmo()
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  //set_cosmological_parameters_to_Joe();
  sprintf(pdeltaparams.runmode,"Halofit");
}
void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  //set_cosmological_parameters_to_Joe();
  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}
void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int cluster_rich_Nbin)
{
  printf("-------------------------------------------\n");
  printf("Initializing Binning\n");
  printf("-------------------------------------------\n");
  
  like.Rmin_bias=Rmin_bias;
  like.Ncl=Ncl;
  like.lmin= lmin; //std=20
  like.lmax= lmax; //15,000
  like.lmax_shear = lmax_shear; //5000
  //compute cluster ell bins acc to 2PCF l-bins
  double ell;
  int i,k=0;
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for(i=0;i<like.Ncl;i++){
    ell=exp(log(like.lmin)+(i+0.5)*logdl);
    if (ell > like.lmax_shear){
      if (k==0) Cluster.l_min = ell;
      k=k+1;
    }
  } 
  Cluster.N200_Nbin = cluster_rich_Nbin;
  Cluster.lbin = k;
  Cluster.l_max = lmax; //clusters go to highly nonlin as std
  printf("%le %le %d\n",Cluster.l_min,Cluster.l_max,Cluster.lbin);
  like.lmax_kappacmb = 2999.;
  
  printf("number of ell bins Ncl: %d\n",like.Ncl);
  printf("minimum ell: %le\n",like.lmin);
  printf("maximum ell: %le\n",like.lmax);
}

void init_binning_real(int Nt, double min, double max)
{
  like.Ntheta= Nt; //std 15
  like.vtmin= min*constants.arcmin; //std 2.5
  like.vtmax= max*constants.arcmin; //std 250.0
}

void init_priors(char *cosmoPrior1, char *cosmoPrior2, char *cosmoPrior3, char *cosmoPrior4)
{
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing priors for marginalization\n");
  printf("---------------------------------------\n");
  
  like.Planck=like.BAO=like.Aubourg_Planck_BAO_SN=like.SN=0;

  if(strcmp(cosmoPrior1,"Planck_BAO_SN_Aubourg")==0)like.Aubourg_Planck_BAO_SN=1;
  if(strcmp(cosmoPrior2,"DES_SN")==0) like.SN=1;
  if(strcmp(cosmoPrior3,"PhotoBAO")==0) like.BAO=1;
  if(strcmp(cosmoPrior4,"Planck")==0) like.Planck=1;  
}


void init_survey(char *surveyname)
{
  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing Survey Parameters\n");
  printf("-------------------------------\n");

  if(strcmp(surveyname,"DES")==0) set_survey_parameters_to_DES();
  if(strcmp(surveyname,"DES_Y1")==0) set_survey_parameters_to_DES_Y1();
  if(strcmp(surveyname,"DES_SV")==0){
    set_survey_parameters_to_DES_SV();
    // sprintf(redshift.shear_REDSHIFT_FILE,"/home/teifler/Dropbox/cosmolike/zdistris/n_of_zs.hist");
    // sprintf(redshift.clustering_REDSHIFT_FILE,"/home/teifler/Dropbox/cosmolike/zdistris/zdistri_redm_SVA");

    sprintf(redshift.shear_REDSHIFT_FILE,"../../zdistris/n_of_zs.hist");
    sprintf(redshift.clustering_REDSHIFT_FILE,"../../zdistris/zdistri_redm_SVA");
    redshift.clustering_photoz=4;
    redshift.shear_photoz = 4;
    redshift.shear_zdistrpar_zmin = 0.01;
    redshift.shear_zdistrpar_zmax = 1.79;
    tomo.shear_Nbin        = 3;
    tomo.shear_Npowerspectra = (int) (tomo.shear_Nbin*(tomo.shear_Nbin+1)/2);
    tomo.shear_zmax[0] = 1.79;
    tomo.shear_zmax[1] = 1.79;
    tomo.shear_zmax[2] = 1.79;
    tomo.shear_zmin[0] = 0.01;
    tomo.shear_zmin[1] = 0.01;
    tomo.shear_zmin[2] = 0.01;
    set_galaxies_DES_SV(0.06);
    set_shear_priors_stage3();
    set_wlphotoz_priors_SV();
  }

  if(strcmp(surveyname,"HSC")==0) set_survey_parameters_to_HSC();  
  if(strcmp(surveyname,"LSST")==0) set_survey_parameters_to_LSST();
  if(strcmp(surveyname,"Euclid")==0) set_survey_parameters_to_Euclid();
  if(strcmp(surveyname,"WFIRST")==0) set_survey_parameters_to_WFIRST();
  printf("Survey set to %s\n",survey.name);
  printf("Survey area: %le deg^2\n",survey.area);
  printf("Source Galaxy Density: %le galaxies/arcmin^2\n",survey.n_gal); 
}


void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *galsample)
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing galaxy samples\n");
  printf("-----------------------------------\n");
  
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",SOURCE_ZFILE);
  printf("PATH TO SOURCE_ZFILE: %s\n",redshift.shear_REDSHIFT_FILE);
  init_source_sample(sourcephotoz);

  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",LENS_ZFILE);
  printf("\n");
  printf("PATH TO LENS_ZFILE: %s\n",redshift.clustering_REDSHIFT_FILE);
  init_lens_sample(lensphotoz,galsample);

  if (strcmp(galsample,"redmagic")==0) set_clphotoz_priors_redmagic();
  if (strcmp(galsample,"benchmark")==0) set_clphotoz_priors_benchmark();
  if (strcmp(galsample,"cmass")==0) set_clphotoz_priors_cmass();
  if (strcmp(galsample,"gold")==0) set_clphotoz_priors_LSST_gold();
  if (strcmp(galsample,"source")==0) set_clphotoz_priors_source();
}

void init_clusters()
{
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing clusters\n");
  printf("-----------------------------------\n");

  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0 || strcmp(survey.name,"HSC")==0) set_clusters_LSST();
  if (strcmp(survey.name,"Euclid")==0 || strcmp(survey.name,"DES")==0)
     set_clusters_DES();

 set_clusterMobs_priors(); 
}


void init_cmb(char * cmbName) {
   printf("\n");
   printf("-----------------------------------\n");
   printf("Initializing CMB\n");
   printf("-----------------------------------\n");
   
   printf("CMB survey: %s\n", cmbName);
   if (strcmp(cmbName, "actpol")==0)
      set_cmb_actpol();
   if (strcmp(cmbName, "advact")==0)
      set_cmb_advact();
   if (strcmp(cmbName, "cmbs4")==0)
      set_cmb_cmbs4();
}

void set_cmb_actpol() {
   sprintf(cmb.name, "actpol");
   cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   cmb.sensitivity = 18.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cov/cmblensrec/actpol/cmblensrecnoise_lmax3000.txt";
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void set_cmb_advact() {
   sprintf(cmb.name, "advact");
   cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   cmb.sensitivity = 10.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cov/cmblensrec/advact/cmblensrecnoise_lmax3000.txt";
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void set_cmb_cmbs4() {
   sprintf(cmb.name, "cmbs4");
   cmb.fwhm = 1. * (constants.pi/180.) / 60.;
   cmb.sensitivity = 1.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cov/cmblensrec/cmbs4/cmblensrecnoise_lmax3000.txt";
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

 
void init_probes(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 

  sprintf(like.probes,"%s",probes);
  if(strcmp(probes,"clusterN")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.clusterN=1;
    printf("Cluster Number Counts computation initialized\n");
  }
  if(strcmp(probes,"clusterN_clusterWL")==0){
    like.Ndata=tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.clusterN=1;
    like.clusterWL=1;
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
  if(strcmp(probes,"all_2pt_clusterN")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
  }

  if(strcmp(probes,"shear_shear")==0){
    like.Ndata=like.Ncl*tomo.shear_Npowerspectra;
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }
  if(strcmp(probes,"pos_pos")==0){
    like.Ndata= like.Ncl*tomo.clustering_Npowerspectra;
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }

  if(strcmp(probes,"all_2pt")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  if(strcmp(probes,"all_2pt_clusterN_clusterWL")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+tomo.cluster_Nbin*Cluster.N200_Nbin+tomo.cgl_Npowerspectra*Cluster.N200_Nbin*Cluster.lbin;
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.clusterN=1;
    like.clusterWL=1;
    printf("%d\n",like.Ndata);
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Cluster Number Counts computation initialized\n");
    printf("Cluster weak lensing computation initialized\n");
  }
   if (strcmp(probes,"LSSxCMB")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: gg, gk, gs, kk, ks, ss\n");
   }
   if (strcmp(probes,"gg_gk_gs")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    printf("Initializing: gg, gk, gs\n");
  }
  if (strcmp(probes,"kk_ks_ss")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Initializing: kk, ks, ss\n");
  }
  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}


void init_probes_real(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 

  sprintf(like.probes,"%s",probes);
  like.shear_shear=0;
  like.shear_pos=0;
  like.pos_pos=0;

  if(strcmp(probes,"shear_shear")==0){
    like.Ndata=like.Ntheta*2*tomo.shear_Npowerspectra;
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }
  
  if(strcmp(probes,"pos_pos")==0){
    like.Ndata= like.Ntheta*tomo.clustering_Npowerspectra;
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }

  if(strcmp(probes,"shear_pos")==0){
    like.Ndata= like.Ntheta*tomo.ggl_Npowerspectra;
    like.shear_pos=1;
    printf("Position-Shear computation initialized\n");
  }

  if(strcmp(probes,"all_2pt")==0){
    like.Ndata=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  
  if(strcmp(probes,"ggl_cl")==0){
    like.Ndata=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  like.Ndata=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);

  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}

void init_data_real(char *COV_FILE, char *MASK_FILE, char *DATA_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.MASK_FILE,"%s",MASK_FILE);
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  init=mask(0,1);

  sprintf(like.COV_FILE,"%s",COV_FILE);
  printf("PATH TO COV: %s\n",like.COV_FILE);
  init=invcov_mask(0,1,1);

  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  init=data_read(0,1);
}

void init_data_inv(char *INV_FILE, char *DATA_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.INV_FILE,"%s",INV_FILE);
  printf("PATH TO INVCOV: %s\n",like.INV_FILE);
  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  init=data_read(0,1);
  init=invcov_read(0,1,1);
}

void init_lens_sample_()
{
  int i,j,n;

  if(strcmp(survey.lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(survey.lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(survey.lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(survey.lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(survey.lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  //  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",survey.lensphotoz,redshift.clustering_photoz);
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    // printf("zmean_lens=%f\n",zmean(i));
    nuisance.bias_zphot_clustering[i]=0.0;
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
    }
  }
  tomo.ggl_Npowerspectra = n;

  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
  // printf("end of lens sample init\n");
}

void init_lens_sample(char *lensphotoz, char *galsample)
{
  if(strcmp(lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  
  if ((redshift.clustering_photoz !=0) && (redshift.clustering_photoz !=1) && (redshift.clustering_photoz !=2) && (redshift.clustering_photoz !=3)) 
  {
    printf("init.c: init_lens_sample: redshift.clustering_photoz = %d not set properly!\nEXIT!\n",redshift.clustering_photoz);
    exit(1);
  }
  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",lensphotoz,redshift.clustering_photoz);
  
  if (strcmp(survey.name,"DES_SV")==0){
    if(strcmp(galsample,"redmagic")==0) {
      init_clphotoz_redmagic();
      set_galaxies_DES_SV(0.06);
    }
  }
  if (strcmp(survey.name,"DES_Y1")==0){ 
    if(strcmp(galsample,"redmagic")==0) {
      init_clphotoz_redmagic();
      set_galaxies_DES_Y1();
    }
    if(strcmp(galsample,"benchmark")==0){ 
      init_clphotoz_benchmark();
      set_galaxies_DES(2.0);
    }
  }
  if (strcmp(survey.name,"DES")==0){ 
    if(strcmp(galsample,"redmagic")==0) {
      init_clphotoz_redmagic();
      set_galaxies_DES(0.15);
    }
    if(strcmp(galsample,"benchmark")==0){ 
      init_clphotoz_benchmark();
      set_galaxies_DES(2.0);
    }
  }

  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0 || strcmp(survey.name,"HSC")==0 || strcmp(survey.name,"Euclid")==0 ){ 
    if(strcmp(galsample,"redmagic")==0){ 
      init_clphotoz_redmagic();   
      set_galaxies_LSST(0.25);
    }
    if(strcmp(galsample,"benchmark")==0){ 
      init_clphotoz_benchmark();
      set_galaxies_LSST(5.0);
    }
    if(strcmp(galsample,"gold")==0){ 
      init_clphotoz_LSST_gold();
      set_galaxies_LSST_gold(40.0);
    }
     if(strcmp(galsample,"cmass")==0){
      init_clphotoz_cmass();
      set_galaxies_CMASS(0.025);
     }
     if(strcmp(galsample,"source")==0){
      init_clphotoz_source();
      set_galaxies_source();
     }
  }
  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
}

void init_source_sample_()
{
  int i;
  if(strcmp(survey.sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(survey.sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(survey.sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(survey.sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(survey.sourcephotoz,"multihisto")==0) redshift.shear_photoz=4;
  //printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",survey.sourcephotoz,redshift.shear_photoz);

  if (redshift.shear_photoz!=4)  set_equal_tomo_bins(tomo.shear_Nbin);
  for (i=0;i<tomo.shear_Nbin; i++)
  {
    //printf("zmean_source=%f\n",zmean_source(i));
    nuisance.bias_zphot_shear[i]=0.0;
  }
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
}


void init_source_sample(char *sourcephotoz)
{
  if(strcmp(sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(sourcephotoz,"multihisto")==0) {
    printf("redshift.shear_photoz=4 not supported\n"); 
    exit(1);
  }
  if ((redshift.shear_photoz !=0) && (redshift.shear_photoz !=1) && (redshift.shear_photoz !=2) && (redshift.shear_photoz !=3)) 
  {
    printf("init.c: init_source_sample: redshift.shear_photoz = %d not set properly!\nEXIT!\n",redshift.shear_photoz);
    exit(1);
  }

  printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",sourcephotoz,redshift.shear_photoz);

  if (strcmp(survey.name,"DES_SV")==0){
    set_equal_tomo_bins(3);
    if ((redshift.shear_photoz==1) || (redshift.shear_photoz==2) || (redshift.shear_photoz==3)){
      init_wlphotoz_stage3();
      set_wlphotoz_priors_stage3();  
    }
  }
  if (strcmp(survey.name,"DES_Y1")==0 || strcmp(survey.name,"DES")==0){
    set_equal_tomo_bins(5);
    if ((redshift.shear_photoz==1) || (redshift.shear_photoz==2) || (redshift.shear_photoz==3)){
      init_wlphotoz_stage3();
      set_wlphotoz_priors_stage3();  
    }
  }
  if (strcmp(survey.name,"Euclid")==0){
    set_equal_tomo_bins(5);
    if ((redshift.shear_photoz==1) || (redshift.shear_photoz==2) || (redshift.shear_photoz==3)){
      init_wlphotoz_stage4();
      set_wlphotoz_priors_stage4();  
    }
  }
  if (strcmp(survey.name,"HSC")==0 ){
    set_equal_tomo_bins(5);
    if ((redshift.shear_photoz==1) || (redshift.shear_photoz==2) || (redshift.shear_photoz==3)){
      init_wlphotoz_stage3();
      set_wlphotoz_priors_stage3();  
    }
  }
  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0) {
    set_equal_tomo_bins(10);
    if ((redshift.shear_photoz==1) || (redshift.shear_photoz==2) || (redshift.shear_photoz==3)){
      init_wlphotoz_stage4();
      set_wlphotoz_priors_stage4();  
    }
  } 
  if (strcmp(survey.name,"LSST")==0 || strcmp(survey.name,"WFIRST")==0 || strcmp(survey.name,"Euclid")==0)  set_shear_priors_stage4();
  if (strcmp(survey.name,"HSC")==0 || strcmp(survey.name,"DES")==0|| strcmp(survey.name,"DES_Y1")==0 || strcmp(survey.name,"DES_SV")==0) set_shear_priors_stage3();
}

void init_baryons()
{
  printf("todo\n");
}

void set_equal_tomo_bins(int Ntomo)
{
  int k,j;
  double frac, zi;
  tomo.shear_Nbin=Ntomo;
  tomo.shear_Npowerspectra=(int) (Ntomo*(Ntomo+1)/2);
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*Ntomo)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

/******init galaxies *******/
//configurations for redmagic as lens galaxies
void set_galaxies_DES_Y1(void)
{
  // Y1 with four lens redshift bins  - redmagic-high_dens for 0..2, redmagic-high_lum for 3
  int i,j,n;
  tomo.clustering_Nbin        = 4;
  tomo.clustering_Npowerspectra = 4;
  tomo.clustering_zmax[0]      = .3;
  tomo.clustering_zmax[1]      = .45;
  tomo.clustering_zmax[2]      = .6;
  tomo.clustering_zmax[3]      = .75;
  
  tomo.clustering_zmin[0]      = 0.15;
  tomo.clustering_zmin[1]      = 0.3;
  tomo.clustering_zmin[2]      = 0.45;
  tomo.clustering_zmin[3]      = 0.6;

  // lens densities estimated from redmagic catalogs
  tomo.n_lens[0] = 0.013;
  tomo.n_lens[1] = 0.034;
  tomo.n_lens[2] = 0.050;
  tomo.n_lens[3] = 0.032;
  

  redshift.clustering_zdistrpar_zmin = 0.1;
  redshift.clustering_zdistrpar_zmax = 0.8;

  printf("\n");
  printf("Lens Sample: DES_Y1 - Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void set_galaxies_DES_SV(double density)
{
  int i,j,n;
  survey.n_lens = density;
  printf("Number density of lens galaxies=%le\n",survey.n_lens);
  tomo.clustering_Nbin        = 2;
  tomo.clustering_Npowerspectra = 2;
  tomo.clustering_zmax[0]      = .6;
  tomo.clustering_zmax[1]      = .6;
  
  tomo.clustering_zmin[0]      = 0.15;
  tomo.clustering_zmin[1]      = 0.15;
  
  redshift.clustering_zdistrpar_zmin = tomo.clustering_zmin[0];
  redshift.clustering_zdistrpar_zmax = tomo.clustering_zmax[tomo.clustering_Nbin-1];

  printf("\n");
  printf("Lens Sample: DES_SV - Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void set_galaxies_DES(double density)
{
  int i,j,n;
  
  redshift.clustering_histogram_zbins=75;

  survey.n_lens = density;
  printf("Number density of lens galaxies=%le\n",survey.n_lens);
  tomo.clustering_Nbin        = 3;
  tomo.clustering_Npowerspectra = 3;

  tomo.clustering_zmax[0]      = .35;
  tomo.clustering_zmax[1]      = .5;
  tomo.clustering_zmax[2]      = .65;
  
  tomo.clustering_zmin[0]      = 0.2;
  tomo.clustering_zmin[1]      = 0.35;
  tomo.clustering_zmin[2]      = 0.5;
  
  redshift.clustering_zdistrpar_zmin = tomo.clustering_zmin[0];
  redshift.clustering_zdistrpar_zmax = tomo.clustering_zmax[tomo.clustering_Nbin-1];

  printf("\n");
  printf("Lens Sample: DES - Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.35+0.15*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
       printf("%d %d %d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  init_BAO_DES();
}
void set_galaxies_LSST_like(double density)
{ //redmagic-like lens sample beyond DES depth
  int i,j,n;
  
  // next 3 lines can perhaps go, beacuse they are being set automatically
  redshift.clustering_zdistrpar_zmin = 0.;
  redshift.clustering_zdistrpar_zmax = 3.5;
  redshift.clustering_histogram_zbins= 75;
  
  survey.n_lens = density; //guestimate of stage 4 lens sample with excellent photo-z
  printf("Number density of lens galaxies=%le\n",survey.n_lens);


  tomo.clustering_Nbin        = 10;
  tomo.clustering_Npowerspectra = 10;
  tomo.clustering_zmax[0] = 0.2;
  tomo.clustering_zmax[1] = 0.4;
  tomo.clustering_zmax[2] = 0.6;
  tomo.clustering_zmax[3] = 0.8;
  tomo.clustering_zmax[4] = 1.0;
  tomo.clustering_zmax[5] = 1.2;
  tomo.clustering_zmax[6] = 1.4;
  tomo.clustering_zmax[7] = 1.6;
  tomo.clustering_zmax[8] = 1.8;
  tomo.clustering_zmax[9] = 2.0;

  tomo.clustering_zmin[0] = 0.05;
  tomo.clustering_zmin[1] = 0.2;
  tomo.clustering_zmin[2] = 0.4;
  tomo.clustering_zmin[3] = 0.6;
  tomo.clustering_zmin[4] = 0.8;
  tomo.clustering_zmin[5] = 1.0;
  tomo.clustering_zmin[6] = 1.2;
  tomo.clustering_zmin[7] = 1.4;
  tomo.clustering_zmin[8] = 1.6;
  tomo.clustering_zmin[9] = 1.8;
  printf("\n");
  printf("Lens Sample: LSST like - Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");

  for (i =0; i < tomo.clustering_Nbin ; i++){
    double z = 0.5*(tomo.clustering_zmin[i] + tomo.clustering_zmax[i]);
    gbias.b[i] = 0.95/growfac(1./(1.+z));
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n; 
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  if (strcmp(survey.name,"LSST")==0) init_BAO_LSST();
  if (strcmp(survey.name,"WFIRST")==0) init_BAO_WFIRST();
} 
void set_galaxies_LSST_gold(double density)
{ //redmagic-like lens sample beyond DES depth
  int i,j,n;
  
  // next 3 lines can perhaps go, beacuse they are being set automatically
  redshift.clustering_zdistrpar_zmin = 0.;
  redshift.clustering_zdistrpar_zmax = 3.5;
  redshift.clustering_histogram_zbins= 75;
  
  survey.n_lens = density; //guestimate of stage 4 lens sample with excellent photo-z
  printf("Number density of lens galaxies=%le\n",survey.n_lens);


  tomo.clustering_Nbin        = 10;
  tomo.clustering_Npowerspectra = 10;
  tomo.clustering_zmax[0] = 0.382940014498;
  tomo.clustering_zmax[1] = 0.522576410401;
  tomo.clustering_zmax[2] = 0.646638132366;
  tomo.clustering_zmax[3] = 0.768971236241;
  tomo.clustering_zmax[4] = 0.897480798253;
  tomo.clustering_zmax[5] = 1.04015552136;
  tomo.clustering_zmax[6] = 1.20895143424;
  tomo.clustering_zmax[7] = 1.42820693242;
  tomo.clustering_zmax[8] = 1.77139821814;
  tomo.clustering_zmax[9] = 3.49999;

  tomo.clustering_zmin[0] = 0.15;
  tomo.clustering_zmin[1] = 0.382940014498;
  tomo.clustering_zmin[2] = 0.522576410401;
  tomo.clustering_zmin[3] = 0.646638132366;
  tomo.clustering_zmin[4] = 0.768971236241;
  tomo.clustering_zmin[5] = 0.897480798253;
  tomo.clustering_zmin[6] = 1.04015552136;
  tomo.clustering_zmin[7] = 1.20895143424;
  tomo.clustering_zmin[8] = 1.42820693242;
  tomo.clustering_zmin[9] = 1.77139821814;

  printf("\n");
  printf("Lens Sample: LSST redmagic- Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");

  for (i =0; i < tomo.clustering_Nbin ; i++){
    double z = 0.5*(tomo.clustering_zmin[i] + tomo.clustering_zmax[i]);
    gbias.b[i] = 0.95/growfac(1./(1.+z));
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n; 
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  if (strcmp(survey.name,"LSST")==0) init_BAO_LSST();
  if (strcmp(survey.name,"WFIRST")==0) init_BAO_WFIRST();
}

void set_galaxies_LSST(double density)
{ //redmagic-like lens sample beyond DES depth
  int i,j,n;
  
  redshift.clustering_zdistrpar_zmin = 0.2;
  redshift.clustering_zdistrpar_zmax = 1.0;
  redshift.clustering_histogram_zbins=75;
  
  survey.n_lens = density; //guestimate of stage 4 lens sample with excellent photo-z
  printf("Number density of lens galaxies=%le\n",survey.n_lens);

  tomo.clustering_Nbin        = 4;
  tomo.clustering_Npowerspectra = 4;
  tomo.clustering_zmax[0]      = .4;
  tomo.clustering_zmax[1]      = .6;
  tomo.clustering_zmax[2]      = .8;
  tomo.clustering_zmax[3]      = 1.;
  
  tomo.clustering_zmin[0]      = 0.2;
  tomo.clustering_zmin[1]      = 0.4;
  tomo.clustering_zmin[2]      = 0.6;
  tomo.clustering_zmin[3]      = 0.8;
  printf("\n");
  printf("Lens Sample: LSST redmagic- Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  printf("\n");
  printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");

  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.35+0.15*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
    }
  }
  tomo.ggl_Npowerspectra = n; 
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
  if (strcmp(survey.name,"LSST")==0) init_BAO_LSST();
  if (strcmp(survey.name,"WFIRST")==0) init_BAO_WFIRST();
}

// !!!!!!!!!!! MANUWARNING: do we want to use HSC as lens sample?

// !!!!!!!!!!! MANUWARNING: set_galaxies_CMASS?
void set_galaxies_CMASS(double density)
{
   int i,j,n;
// MANUWARNING: these parameters?
   redshift.clustering_zdistrpar_zmin = 0.4; //0.2;
   redshift.clustering_zdistrpar_zmax = 0.7; //1.0;
   redshift.clustering_histogram_zbins = 30;
   
   survey.n_lens = density; //guestimate of stage 4 lens sample with excellent photo-z
   printf("Number density of source galaxies=%le\n",survey.n_lens);
   
   tomo.clustering_Nbin        = 2;
   tomo.clustering_Npowerspectra = 2;
   tomo.clustering_zmax[0]      = 0.57;
   tomo.clustering_zmax[1]      = 0.7;
   
   tomo.clustering_zmin[0]      = 0.4;
   tomo.clustering_zmin[1]      = 0.57;
   printf("\n");
   printf("Lens Sample: CMASS - Tomographic Bin limits:\n");
   for (i =0; i<tomo.clustering_Nbin ; i++){
      printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
   }
   printf("\n");
   printf("Setting Galaxy Bias - Passive Evolution in z-bins:\n");
   
   for (i =0; i < tomo.clustering_Nbin ; i++){
      gbias.b[i] = 2.;
      printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
   }
   n = 0;
   for (i = 0; i<tomo.clustering_Nbin; i++){
      for(j = 0; j<tomo.shear_Nbin;j++){
         n += test_zoverlap(i,j);
      }
   }
   tomo.ggl_Npowerspectra = n;
   printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}


void set_galaxies_source(void)
{
  // Y1 with four lens redshift bins  - redmagic-high_dens for 0..2, redmagic-high_lum for 3
  int i,j,n;
  tomo.clustering_Nbin        = tomo.shear_Nbin;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    tomo.clustering_zmax[i]      = tomo.shear_zmax[i];
    tomo.clustering_zmin[i]      = tomo.shear_zmin[i];
  }
  tomo.clustering_zmin[0]=0.2; //multi-probe NG covs are very likely ill-conditioned if lenses at very low redshift is included 
  survey.n_lens = survey.n_gal; 

  redshift.clustering_zdistrpar_zmin = redshift.shear_zdistrpar_zmin;
  redshift.clustering_zdistrpar_zmax = redshift.shear_zdistrpar_zmax;

  printf("\n");
  printf("Lens Sample: Source - Tomographic Bin limits:\n");
  for (i =0; i < tomo.clustering_Nbin ; i++){
    printf("min=%le max=%le\n",tomo.clustering_zmin[i],tomo.clustering_zmax[i]);
  }
  gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}






/*********** set cluster parameters for DES forecasts ********/
void set_clusters_LSST(){
  int i,j;
  //N200-M relationship from Rykoff et al. 2012 (http://iopscience.iop.org/0004-637X/746/2/178/pdf/0004-637X_746_2_178.pdf) - Eq. B4 for \Delta = 200 mean version
  nuisance.cluster_Mobs_lgM0 = 1.72+log(1.e+14*0.7); //back to Msun/h instead of Msun/h70 normalization
  nuisance.cluster_Mobs_sigma = 0.25;
  nuisance.cluster_Mobs_alpha = 1.08;
  nuisance.cluster_Mobs_beta = 0.0;
  nuisance.cluster_Mobs_N_pivot = 60.;
  // completeness - adjusts later
  nuisance.cluster_completeness[0] = 0.9;
  nuisance.cluster_completeness[1] = 0.9;
  nuisance.cluster_completeness[2] = 0.9;
  nuisance.cluster_completeness[3] = 0.9;
  //no miscentering so far
  nuisance.cluster_centering_f0 = 1.0;
  nuisance.cluster_centering_alpha = 0;
  nuisance.cluster_centering_sigma = 0;
  nuisance.cluster_centering_M_pivot = 1.e+14;
  printf("%e %e %e %e\n",nuisance.cluster_Mobs_lgM0,nuisance.cluster_Mobs_alpha,nuisance.cluster_Mobs_beta,nuisance.cluster_Mobs_N_pivot);
  tomo.cluster_Nbin = 4; // number of cluster redshift bins
  tomo.cluster_zmin[0] = 0.2;
  tomo.cluster_zmax[0] = 0.4;
  tomo.cluster_zmin[1] = 0.4;
  tomo.cluster_zmax[1] = 0.6;
  tomo.cluster_zmin[2] = 0.6;
  tomo.cluster_zmax[2] = 0.8;
  tomo.cluster_zmin[3] = 0.8;
  tomo.cluster_zmax[3] = 1.0;
  tomo.cgl_Npowerspectra = 0;// number of cluster-lensing tomography combinations
  for (i = 0; i < tomo.cluster_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      tomo.cgl_Npowerspectra += test_zoverlap_c(i,j);
    }
  }

  Cluster.N200_min = 10.;
  Cluster.N200_max = 220.;
  Cluster.N200_Nbin = 7;
  //upper bin boundaries - note that bin boundaries need to be integers!
  int Nlist[8] = {14,20,30,45,70,120,Cluster.N200_max};
  Cluster.N_min[0] = Cluster.N200_min;
  Cluster.N_max[0] = Nlist[0];
  for (i = 1; i < Cluster.N200_Nbin; i++){
    Cluster.N_min[i] = Nlist[i-1];
    Cluster.N_max[i] = Nlist[i];
  }
 for (i = 0; i < Cluster.N200_Nbin; i++){
    printf ("Richness bin %d: %e - %e (%e Msun/h - %e Msun/h), N(z = 0.3) = %e, N(z = 0.7) = %e\n", i,Cluster.N_min[i],Cluster.N_max[i],exp(lgM_obs(Cluster.N_min[i], 0.75)),exp(lgM_obs(Cluster.N_max[i], 0.75)),N_N200(0,i),N_N200(2,i));
  }
  printf("Clusters set to LSST\n");
  printf("Clusters cgl_Npowerspectra=%d\n",tomo.cgl_Npowerspectra);
}

void set_clusters_DES(){
  int i,j;
  //N200-M relationship from Rykoff et al. 2012 (http://iopscience.iop.org/0004-637X/746/2/178/pdf/0004-637X_746_2_178.pdf) - Eq. B4 for \Delta = 200 mean version
  nuisance.cluster_Mobs_lgM0 = 1.72+log(1.e+14*0.7); //back to Msun/h instead of Msun/h70 normalization
  nuisance.cluster_Mobs_sigma = 0.25;
  nuisance.cluster_Mobs_alpha = 1.08;
  nuisance.cluster_Mobs_beta = 0.0;
  nuisance.cluster_Mobs_N_pivot = 60.;
  // completeness - adjusts later
  nuisance.cluster_completeness[0] = 0.9;
  nuisance.cluster_completeness[1] = 0.9;
  nuisance.cluster_completeness[2] = 0.9;
  nuisance.cluster_completeness[3] = 0.9;
  //no miscentering so far
  nuisance.cluster_centering_f0 = 1.0;
  nuisance.cluster_centering_alpha = 0;
  nuisance.cluster_centering_sigma = 0;
  nuisance.cluster_centering_M_pivot = 1.e+14;
  printf("%e %e %e %e\n",nuisance.cluster_Mobs_lgM0,nuisance.cluster_Mobs_alpha,nuisance.cluster_Mobs_beta,nuisance.cluster_Mobs_N_pivot);
  tomo.cluster_Nbin = 3; // number of cluster redshift bins
  tomo.cluster_zmin[0] = 0.2;
  tomo.cluster_zmax[0] = 0.4;
  tomo.cluster_zmin[1] = 0.4;
  tomo.cluster_zmax[1] = 0.6;
  tomo.cluster_zmin[2] = 0.6;
  tomo.cluster_zmax[2] = 0.8;
  tomo.cgl_Npowerspectra = 0;// number of cluster-lensing tomography combinations
  for (i = 0; i < tomo.cluster_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      tomo.cgl_Npowerspectra += test_zoverlap_c(i,j);
    }
  }
  printf("Cluster Lensing Power Spectra %d\n",tomo.cgl_Npowerspectra);
  Cluster.N200_min = 10.;
  Cluster.N200_max = 220.;
  Cluster.N200_Nbin = 7;
  //upper bin boundaries - note that bin boundaries need to be integers!
  int Nlist[8] = {14,20,30,45,70,120,Cluster.N200_max};
  Cluster.N_min[0] = Cluster.N200_min;
  Cluster.N_max[0] = Nlist[0];
  for (i = 1; i < Cluster.N200_Nbin; i++){
    Cluster.N_min[i] = Nlist[i-1];
    Cluster.N_max[i] = Nlist[i];
  }
  //  for (i = 0; i < Cluster.N200_Nbin; i++){
  //    printf ("Richness bin %d: %e - %e (%e Msun/h - %e Msun/h), N(z = 0.3) = %e, N(z = 0.7)\n", i,Cluster.N_min[i],Cluster.N_max[i],exp(lgM_obs(Cluster.N_min[i], 0.75)),exp(lgM_obs(Cluster.N_max[i], 0.75)),N_N200(0,i),N_N200(2,i));
  //  }
  printf("Clusters set to DES\n");
}

void set_galaxies_benchmark_LSST()
{ 
  printf("todo\n");
}

void set_galaxies_benchmark_DES()
{ 
  printf("todo\n");
}

void init_Pdelta(char *model,double nexp,double A_factor)
{  
  sprintf(pdeltaparams.runmode,"%s",model);
  pdeltaparams.DIFF_n=nexp;
  pdeltaparams.DIFF_A=A_factor;
}


void init_IA(char *model,char *lumfct)
{  
  if(strcmp(lumfct,"GAMA")==0) set_LF_GAMA();
  else if(strcmp(lumfct,"DEEP2")==0) set_LF_DEEP2();
  else {
    printf("init.c:init_IA: %s lumfct not defined\n",lumfct);
    printf("USING GAMA LF INSTEAD\n");
    set_LF_GAMA();
  }
  printf("SET LUMINOSITY FUNCTION=%s\n",lumfct);
  
  nuisance.oneplusz0_ia=1.3; 
  //z0=0.3 is arbitrary pivot redshift J11 p18
  nuisance.c1rhocrit_ia=0.0134; 
  // J11 p.8
  
  if(strcmp(model,"none")==0)  like.IA=0;
  else if(strcmp(model,"NLA_HF")==0)  like.IA=1;
  else if(strcmp(model,"lin")==0)  like.IA=2;
  else{
    printf("init.c:init_IA: %s IA model not defined\n",model);
    exit(1);
  }
  printf("SET IA MODEL=%s\n",model);
  set_ia_priors();
  log_like_f_red();
}


void init_wlphotoz_stage3()
{
  int i;
  printf("\n");
  printf("Source sample: stage 3 photoz uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.08; //acc Troxel
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
  } 
}


void init_wlphotoz_stage4()
{
  int i;
  printf("\n");
  printf("Source sample: stage 4 photoz uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=0.05; 
    printf("nuisance.bias_zphot_shear[%d]=%le\n",i,nuisance.bias_zphot_shear[i]);
    printf("nuisance.sigma_zphot_shear[%d]=%le\n",i,nuisance.sigma_zphot_shear[i]);
  }
}


void init_clphotoz_redmagic()
{
  int i;
  printf("\n");
  printf("Galaxy sample redmagic photoz uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.01; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
  }  
}

void init_clphotoz_LSST_gold()
{
  int i;
  printf("\n");
  printf("Galaxy sample LSST Gold photoz uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.03; 
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]);
  }  
}

void init_clphotoz_benchmark()
{
  int i;
  printf("\n");
  printf("Galaxy sample benchmark photoz uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.04;
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]); 
  }
}

void init_clphotoz_cmass()
{
  int i;
  printf("\n");
  printf("Galaxy sample cmass redshift uncertainty initialized\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=0.001;
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]); 
  }
}

void init_clphotoz_source()
{
  int i;
  printf("\n");
  printf("Lens sample initialized with same parameters as source sample\n");
  for (i=0;i<10; i++){
    nuisance.bias_zphot_clustering[i]=nuisance.bias_zphot_shear[i];
    nuisance.sigma_zphot_clustering[i]=nuisance.sigma_zphot_shear[i];
    printf("nuisance.bias_zphot_clustering[%d]=%le\n",i,nuisance.bias_zphot_clustering[i]);
    printf("nuisance.sigma_zphot_clustering[%d]=%le\n",i,nuisance.sigma_zphot_clustering[i]); 
  }
}

void init_HOD_rm(){
  set_HOD_redmagic_priors();
  like.Rmin_bias = 0.1;//use halo+HOD model down to 100 kpc/h
  redm.parameterization = 0; //Zehavi et al. 2011 HOD parameterization
  redm.cg = 1.0;
  redm.fc = 0.2;
  redm.hod[0] = 12.1;
  redm.hod[1] = 0.4;
  redm.hod[2] = 13.65;
  redm.hod[3] = 12.2;
  redm.hod[4] = 1.0;
}

