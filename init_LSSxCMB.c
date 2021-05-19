double invcov_read(int READ, int ci, int cj);
double data_read(int READ, int ci);
double bary_read(int READ, int PC, int cj);
void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE);
void init_priors(double M_Prior, double SigZ_source, double DeltaZ_source_Prior, double SigZ_source_Prior, double SigZ_lens, double DeltaZ_lens_Prior, double SigZ_lens_Prior, double A_IA_Prior, double betaIA_Prior, double etaIA_Prior, double etaZIA_Prior, double Q1_Prior, double Q2_Prior, double Q3_Prior);
void init_survey(char *surveyname, double nsource, double nlens, double area);
void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *tomo_binning_source, char *tomo_binning_lens);
void init_cosmo_runmode(char *runmode);
void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens);
void init_probes(char *probes);


void init_lens_sample(char *lensphotoz, char *tomo_binning_lens);
void init_source_sample(char *sourcephotoz, char *tomo_binning_source);

void set_galaxies_source();
void set_lens_galaxies_LSSTgoldsample();
void set_galaxies_WFIRST_SN10();
void set_galaxies_lens_as_source();

void set_equal_tomo_bins();
void init_IA(char *model,char *lumfct);

void init_cmb(char * cmbName);
void set_cmb_cmbs4();
void set_cmb_so_Y5();
void set_cmb_so_Y1();

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

  p = line;
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


void init_data_inv_bary(char *INV_FILE, char *DATA_FILE, char *BARY_FILE)
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
  sprintf(like.BARY_FILE,"%s",BARY_FILE);
  printf("PATH TO BARYONS: %s\n",like.BARY_FILE);
  init=data_read(0,1);
  init=bary_read(0,1,1);
  init=invcov_read(0,1,1);
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

double bary_read(int READ, int PC, int cj)
{
  int i,j,intspace, N_PC=6;
  static double **bary =0;

  if(READ==0 || bary == 0){
    bary   = create_double_matrix(0, N_PC-1, 0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.BARY_FILE,"r");
    for (i=0;i<like.Ndata; i++){
      fscanf(F,"%le %le %le %le %le %le\n",&bary[0][i],&bary[1][i],&bary[2][i],&bary[3][i],&bary[4][i],&bary[5][i]);  
    }
    fclose(F);
    printf("FINISHED READING BARYON MATRIX\n");
  }    
  return bary[PC][cj];
}


void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_binning_fourier(int Ncl, double lmin, double lmax, double lmax_shear, double Rmin_bias, int Ntomo_source, int Ntomo_lens)
{
  printf("-------------------------------------------\n");
  printf("Initializing Binning\n");
  printf("-------------------------------------------\n");
  
  like.Rmin_bias=Rmin_bias;
  like.Ncl=Ncl;
  like.lmin= lmin; //std=20
  like.lmax= lmax; //15,000
  like.lmax_shear = lmax_shear; //5000
  tomo.shear_Nbin=Ntomo_source;
  tomo.clustering_Nbin=Ntomo_lens;
  double ell;
  int i,k=0;
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  for(i=0;i<like.Ncl;i++){
    ell=exp(log(like.lmin)+(i+0.5)*logdl);
  } 
  
  printf("number of ell bins Ncl: %d\n",like.Ncl);
  printf("minimum ell: %le\n",like.lmin);
  printf("maximum ell: %le\n",like.lmax);
}


void init_priors(double M_Prior, double SigZ_source, double DeltaZ_source_Prior, double SigZ_source_Prior, double SigZ_lens, double DeltaZ_lens_Prior, double SigZ_lens_Prior, double A_ia_Prior, double beta_ia_Prior, double eta_ia_Prior, double etaZ_ia_Prior, double Q1_Prior, double Q2_Prior, double Q3_Prior)
{
  int i;

  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing observational priors for marginalization\n");
  printf("---------------------------------------\n");
  
  for (i=0;i<tomo.shear_Nbin; i++){
    prior.shear_calibration_m[i][0] = 0.0;
    prior.shear_calibration_m[i][1] = M_Prior;
  }
  like.shearcalib=1;
  
  for (i=0;i<tomo.shear_Nbin; i++){
    nuisance.bias_zphot_shear[i]=0.0;
    nuisance.sigma_zphot_shear[i]=SigZ_source; 
    // center of Gaussian priors
    prior.bias_zphot_shear[i][0]=nuisance.bias_zphot_shear[i];
    prior.sigma_zphot_shear[i][0]=nuisance.sigma_zphot_shear[i]; 
    // rms width of Gaussian priors
    prior.bias_zphot_shear[i][1] = DeltaZ_source_Prior;
    prior.sigma_zphot_shear[i][1]= SigZ_source_Prior;
  }
  like.wlphotoz=1;

  for (i=0;i<tomo.clustering_Nbin; i++){
    nuisance.bias_zphot_clustering[i]=0.0;
    nuisance.sigma_zphot_clustering[i]=SigZ_lens; 
    // center of Gaussian priors
    prior.bias_zphot_clustering[i][0]=nuisance.bias_zphot_clustering[i];
    prior.sigma_zphot_clustering[i][0]=nuisance.sigma_zphot_clustering[i];
    // rms width of Gaussian priors
    prior.bias_zphot_clustering[i][1] = DeltaZ_lens_Prior;
    prior.sigma_zphot_clustering[i][1]= SigZ_lens_Prior;
  }
  like.clphotoz=1;

  prior.A_ia[0]=5.92; 
  prior.A_ia[1]=A_ia_Prior; 
  
  prior.beta_ia[0]=1.1; 
  prior.beta_ia[1]=beta_ia_Prior; 
  
  prior.eta_ia[0]=-0.47; 
  prior.eta_ia[1]=eta_ia_Prior; 
  
  prior.eta_ia_highz[0]=0.0; 
  prior.eta_ia_highz[1]=etaZ_ia_Prior; 
  like.IA=1;

  prior.bary_Q1[0]=0.0; 
  prior.bary_Q1[1]=Q1_Prior; 

  prior.bary_Q2[0]=0.0; 
  prior.bary_Q2[1]=Q2_Prior; 

  prior.bary_Q3[0]=0.0; 
  prior.bary_Q3[1]=Q3_Prior; 
  like.baryons=1;

#ifdef NOMPP
  printf("\n");
  printf("---------------------------------------\n");
  printf("Shear Calibration Prior\n");
  printf("Mean=%le, Sigma=%le\n",prior.shear_calibration_m[0][0],prior.shear_calibration_m[i][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Weak Lensing\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_shear[0][0],prior.bias_zphot_shear[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_shear[0][0],prior.sigma_zphot_shear[0][1]); 

  printf("\n");
  printf("---------------------------------------\n");
  printf("Photo-z priors Clustering\n");
  printf("Delta_z=%le, Sigma (Delta_z)=%le\n",prior.bias_zphot_clustering[0][0],prior.bias_zphot_clustering[0][1]);
  printf("Sigma_z=%le, Sigma (Sigma_z)=%le\n",prior.sigma_zphot_clustering[0][0],prior.sigma_zphot_clustering[0][1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("IA Priors\n");
  printf("A_IA=%le, A_IA_Prior=%le\n",prior.A_ia[0],prior.A_ia[1]);
  printf("beta_ia=%le, betaIA_Prior=%le\n",prior.beta_ia[0],prior.beta_ia[1]);
  printf("eta_ia=%le, etaIA_Prior=%le\n",prior.eta_ia[0],prior.eta_ia[1]);
  printf("eta_ia_highz=%le, etaZIA_Prior=%le\n",prior.eta_ia_highz[0],prior.eta_ia_highz[1]);

  printf("\n");
  printf("---------------------------------------\n");
  printf("Baryon Priors\n");
  printf("Q1=%le, Sigma (Q1)=%le\n",prior.bary_Q1[0],prior.bary_Q1[1]);
  printf("Q2=%le, Sigma (Q2)=%le\n",prior.bary_Q2[0],prior.bary_Q2[1]);
  printf("Q3=%le, Sigma (Q3)=%le\n",prior.bary_Q3[0],prior.bary_Q3[1]);
#endif
}
 

void init_survey(char *surveyname, double nsource, double nlens, double area)
{
#ifdef NOMPP
  printf("\n");
  printf("-------------------------------\n");
  printf("Initializing Survey Parameters\n");
  printf("-------------------------------\n");
#endif

  survey.area   = area;
  survey.n_gal   = nsource;
  survey.n_lens=nlens;
  survey.sigma_e   = 0.37;
  sprintf(survey.name,"%s",surveyname);
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
  survey.m_lim=24.5;
#ifdef NOMPP
  printf("Survey set to %s\n",survey.name);
  printf("Survey area: %le deg^2\n",survey.area);
  printf("Source Galaxy Density: %le galaxies/arcmin^2\n",survey.n_gal);
#endif
}


void init_galaxies(char *SOURCE_ZFILE, char *LENS_ZFILE, char *lensphotoz, char *sourcephotoz, char *tomo_binning_source, char *tomo_binning_lens)
{
#ifdef NOMPP
  printf("\n");
  printf("-----------------------------------\n");
  printf("Initializing galaxy samples\n");
  printf("-----------------------------------\n");
#endif
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",SOURCE_ZFILE);
#ifdef NOMPP
  printf("PATH TO SOURCE_ZFILE: %s\n",redshift.shear_REDSHIFT_FILE);
#endif
  init_source_sample(sourcephotoz,tomo_binning_source);
  
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",LENS_ZFILE);
#ifdef NOMPP
  printf("\n");
  printf("PATH TO LENS_ZFILE: %s\n",redshift.clustering_REDSHIFT_FILE);
#endif
  init_lens_sample(lensphotoz,tomo_binning_lens);
}



void init_probes(char *probes)
{
#ifdef NOMPP
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing Probes\n");
  printf("------------------------------\n"); 
  printf("like.Ncl=%d\n",like.Ncl);
  printf("tomo.shear_Npowerspectra=%d\n",tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra=%d\n",tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra=%d\n",tomo.clustering_Npowerspectra);
#endif

  sprintf(like.probes,"%s",probes);
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
  if(strcmp(probes,"ggl_cl")==0){
    like.Ndata=like.Ncl*(tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  if(strcmp(probes,"3x2pt")==0){
    like.Ndata=like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  if(strcmp(probes,"5x2pt")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-Position computation initialized\n");
  }
  if(strcmp(probes,"6x2pt")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra+1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-Position computation initialized\n");
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"gg_gk_gs")==0) {
    like.Ndata = like.Ncl * (2*tomo.clustering_Nbin+tomo.ggl_Npowerspectra);
    like.pos_pos = 1;
    like.gk = 1;
    like.shear_pos = 1;
    printf("Position-Position computation initialized\n");
    printf("Position-Shear computation initialized\n");
    printf("Position-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"kk_ks_ss")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.shear_Npowerspectra);
    like.kk = 1;
    like.ks = 1;
    like.shear_shear = 1;
    printf("Shear-Shear computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"gk_ks_kk")==0) {
    like.Ndata = like.Ncl * (1+tomo.shear_Nbin+tomo.clustering_Nbin);
    like.kk = 1;
    like.ks = 1;
    like.gk = 1;
    printf("Position-CMBkappa computation initialized\n");
    printf("CMBkappa-Shear computation initialized\n");
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  if(strcmp(probes,"kk")==0) {
    like.Ndata = like.Ncl;
    like.kk = 1;
    printf("CMBkappa-CMBkappa computation initialized\n");
  }
  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}



void init_lens_sample(char *lensphotoz, char *tomo_binning_lens)
{
  if(strcmp(lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  if(strcmp(lensphotoz,"outlier_sim")==0) redshift.clustering_photoz=5;
  if(strcmp(lensphotoz,"outlier_model")==0) redshift.clustering_photoz=6;
  
  if ((redshift.clustering_photoz !=0) && (redshift.clustering_photoz !=1) && (redshift.clustering_photoz !=2) && (redshift.clustering_photoz !=3) && (redshift.clustering_photoz !=5) && (redshift.clustering_photoz !=6)) 
  {
    printf("init.c: init_lens_sample: redshift.clustering_photoz = %d not set properly!\nEXIT!\n",redshift.clustering_photoz);
    exit(1);
  }
  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",lensphotoz,redshift.clustering_photoz);
  
  if(strcmp(tomo_binning_lens,"WF_SN10")==0){
    set_galaxies_WFIRST_SN10();
  }
  if(strcmp(tomo_binning_lens,"LSST_gold")==0){
    set_lens_galaxies_LSSTgoldsample();
  }
  if(strcmp(tomo_binning_lens,"lens=src")==0){
    set_galaxies_lens_as_source();
  }
  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
}



void init_source_sample(char *sourcephotoz, char *tomo_binning_source)
{
  if(strcmp(sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(sourcephotoz,"multihisto")==0) {
    printf("redshift.shear_photoz=4 not supported\n"); 
    exit(1);
  }
  if(strcmp(sourcephotoz,"outlier_sim")==0) redshift.shear_photoz=5;
  if(strcmp(sourcephotoz,"outlier_model")==0) redshift.shear_photoz=6;
  if ((redshift.shear_photoz !=0) && (redshift.shear_photoz !=1) && (redshift.shear_photoz !=2) && (redshift.shear_photoz !=3) && (redshift.shear_photoz !=5) && (redshift.shear_photoz !=6)) 
  {
    printf("init.c: init_source_sample: redshift.shear_photoz = %d not set properly!\nEXIT!\n",redshift.shear_photoz);
    exit(1);
  }

  printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",sourcephotoz,redshift.shear_photoz);
  if(strcmp(tomo_binning_source,"source_std")==0)set_galaxies_source();
}


void set_galaxies_source()
{
  int k,j;
  double frac, zi;
  
  tomo.shear_Npowerspectra=(int) (tomo.shear_Nbin*(tomo.shear_Nbin+1)/2);
  
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;printf("redshift.shear_zdistrpar_zmax,min: %le, %le\n", redshift.shear_zdistrpar_zmax, redshift.shear_zdistrpar_zmin);
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
#ifdef NOMPP
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
#endif

  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.shear_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
#ifdef NOMPP
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
#endif
  }
#ifdef NOMPP
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
#endif
  free_double_vector(sum,0,zbins);
}

void set_galaxies_lens_as_source()
{
  int i,k,j,n;
  double frac, zi;
  
  tomo.clustering_Npowerspectra=tomo.shear_Nbin;
  tomo.clustering_Nbin = tomo.shear_Nbin;
  redshift.clustering_zdistrpar_zmin = redshift.shear_zdistrpar_zmin;
  redshift.clustering_zdistrpar_zmax = redshift.shear_zdistrpar_zmax;

  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.clustering_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = redshift.shear_zdistrpar_zmax;
#ifdef NOMPP
  printf("\n");
  printf("Lens as Source - Tomographic Bin limits:\n");
#endif

  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.clustering_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
#ifdef NOMPP
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
#endif
  }
#ifdef NOMPP
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
#endif
  free_double_vector(sum,0,zbins);

  gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 0.95/(growfac(1./(1.+(tomo.clustering_zmax[i]+tomo.clustering_zmin[i]/2.)))/growfac(1.));
    // gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n=0;

  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void set_lens_galaxies_LSSTgoldsample()
{
  int i,j,n,k;
  double frac, zi;
  redshift.clustering_zdistrpar_zmin = 0.01;
  redshift.clustering_zdistrpar_zmax = 1.5;
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  tomo.clustering_zmin[0] = 0.2;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = 1.2;

  int zbins =2000;
  double da = (tomo.clustering_zmax[tomo.clustering_Nbin-1]-tomo.clustering_zmin[0])/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = tomo.clustering_zmin[0]; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+pf_histo(zi, NULL);
  }
  printf("\n");
  printf("Lens Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = tomo.clustering_zmin[0]+j*da;
    tomo.clustering_zmin[k+1] = tomo.clustering_zmin[0]+j*da;
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
    gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin ; i++){
    gbias.b[i] = 0.95/(growfac(1./(1.+(tomo.clustering_zmax[i]+tomo.clustering_zmin[i]/2.)))/growfac(1.));
    //gbias.b[i] = 1.3+0.1*i;
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
  }
  n=0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
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


/************ CMB Settings ***********/
void init_cmb(char * cmbName) {
   printf("\n");
   printf("-----------------------------------\n");
   printf("Initializing CMB\n");
   printf("-----------------------------------\n");
   
   printf("CMB survey: %s\n", cmbName);
   if (strcmp(cmbName, "cmbs4")==0)
      set_cmb_cmbs4();
   if (strcmp(cmbName, "so_Y1")==0)
      set_cmb_so_Y1();
   if (strcmp(cmbName, "so_Y5")==0)
      set_cmb_so_Y5();
}


void set_cmb_cmbs4() {
   sprintf(cmb.name, "cmbs4");
   cmb.fwhm = 1. * (constants.pi/180.) / 60.;
   cmb.sensitivity = 1.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cmblensrec/cmbs4/cmb_lmax3000.txt";
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

void set_cmb_so_Y5() {
   sprintf(cmb.name, "so_Y5");
   // cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   // cmb.sensitivity = 18.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cmblensrec/so/YEAR5_2colformat_nlkk_v3_1_0deproj0_SENS1_fsky0p4_it_lT30-3000_lP30-5000.dat";
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}
void set_cmb_so_Y1() {
   sprintf(cmb.name, "so_Y1");
   // cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
   // cmb.sensitivity = 18.*(constants.pi/180.)/60.;
   cmb.pathLensRecNoise = "./cmblensrec/so/YEAR1_nlkk_SOlike_y1_tt_SENS1_qe_fsky0p4_lT30-3000.dat";
   like.lmax_kappacmb = 2999.;
   printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
}

// void set_cmb_so_gold() {
//    sprintf(cmb.name, "so_gold");
//    // cmb.fwhm = 1.4 * (constants.pi/180.) / 60.;
//    // cmb.sensitivity = 18.*(constants.pi/180.)/60.;
//    cmb.pathLensRecNoise = "./cmblensrec/so/so_gold_nlkk_lmax3000.txt";
//    like.lmax_kappacmb = 2999.;
//    printf("path for CMB lens noise: %s\n", cmb.pathLensRecNoise);
// }



////// For WFIRSTxSO forecast
void set_galaxies_WFIRST_SN10()
{
  int k,j,n,i;
  double frac, zi;
  redshift.clustering_zdistrpar_zmin = 0.25;
  redshift.clustering_zdistrpar_zmax = 4.0;
  tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
  tomo.clustering_zmin[0] = redshift.clustering_zdistrpar_zmin;
  tomo.clustering_zmax[tomo.clustering_Nbin-1] = redshift.clustering_zdistrpar_zmax;
  pf_histo(zi, NULL);
  int zbins =2000;
  double da = (redshift.clustering_zdistrpar_zmax-redshift.clustering_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.clustering_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+pf_histo(zi, NULL);
  }

#ifdef NOMPP
  printf("\n");
  printf("Lens Sample - Tomographic Bin limits:\n");
#endif
  for(k=0;k<tomo.clustering_Nbin-1;k++){
    frac=(k+1.)/(1.*tomo.clustering_Nbin)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.clustering_zmax[k] = redshift.clustering_zdistrpar_zmin+j*da;
    tomo.clustering_zmin[k+1] = redshift.clustering_zdistrpar_zmin+j*da;
#ifdef NOMPP
    printf("min=%le max=%le\n",tomo.clustering_zmin[k],tomo.clustering_zmax[k]);
#endif
  }
#ifdef NOMPP
  printf("min=%le max=%le\n",tomo.clustering_zmin[tomo.clustering_Nbin-1],tomo.clustering_zmax[tomo.clustering_Nbin-1]);
  printf("redshift.clustering_zdistrpar_zmin=%le max=%le\n",redshift.clustering_zdistrpar_zmin,redshift.clustering_zdistrpar_zmax);
#endif
  free_double_vector(sum,0,zbins);
    gbias.b1_function = &b1_per_bin;
  for (i =0; i < tomo.clustering_Nbin; i++){
    gbias.b[i] = 1.3+0.1*i;
#ifdef NOMPP
    printf("Bin %d: galaxy bias=%le\n",i,gbias.b[i]);
#endif
  }
  n=0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
#ifdef NOMPP
      printf("GGL combinations zl=%d zs=%d accept=%d\n",i,j,test_zoverlap(i,j));
#endif
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}