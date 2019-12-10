double cov_NG_shearxkappaCMB_shearxkappaCMB (double l1,double l2, int z1, int z2);//super-sample + halo trispectrum covariance between powerspectra kappa(z1)-kappa(CMB) and kappa(z2)-kappa(CMB)
double cov_G_shearxkappaCMB_shearxkappaCMB (double l, double delta_l, int z1, int z2);//Gaussian covariance between powerspectra kappa(z1)-kappa(CMB) and kappa(z2)-kappa(CMB)

/********** utility routines ****************/

// Reads in the noise N_mv for mv quadratic estimator for d,
// without reduction due to average over ell-bin.
// Returns value at nearest upper ell
// MANUWARNING: switch between expt, memory leak, scientific format
double kappa_reconstruction_noise(double l){
// !!!!!!!!!!!!!!!!!! how to switch between different surveys? !!!!!!!
// !!!!!!!!!!!!!!!!!! absolute path? !!!!!!!!!!!!!!!
//   char * path = "../top-level/cmb/cov/cmblensrec/plancksmica/cmblensrecnoise_lmax3000.txt";
   char * path = "../top-level/cmb/cov/cmblensrec/actpol/cmblensrecnoise_lmax3000.txt";
//   char * path = "../top-level/cmb/cov/cmblensrec/advact/cmblensrecnoise_lmax3000.txt";
//   char * path = "../top-level/cmb/cov/cmblensrec/cmbs4/cmblensrecnoise_lmax3000.txt";
   
   // if first time
   if (data==0) {
      // count lines
      int nEll = line_count(path);
      // nb of pairs {TT,TT}, {TT,TE}, etc
      int nObs = 16;
      
      // declare arrays
      double **cov[nObs][nEll];
      // allocate ell and Nmv
// !!!!!!!!!!!! these arrays are never freed... memory leak!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double *ell=create_double_vector(0, nEll-1);
      double *noise=create_double_vector(0, nEll-1);
      
      // read each line
      FILE *file = fopen("Integers.txt", "r");
      int iEll;
      for (iEll=0; iEll<nEll; iEll++) {
         // read ell
// !!!!!!!!!!!! is %le the right format? My nbs have 19 sig figs, then exponent !!!!!!!!
         fscanf(file, "%le", &ell[iEll]);
         // read covariances
         int iObs;
         for (iObs=1; iObs<nObs; iObs++) {
            fscanf(file, "%le", &cov[iObs][iEll])
         }
         // keep N_mv
         noise[iEll] = cov[nObs-1][iEll];
      }
      fclose(file);
   }
   
   // if l is in the range
   if ((l>=ell[0]) &&(l<=ell[nEll-1])){
      // find value of ell just above l
      int iEll = 0;
      while (ell[iEll] < l) {
         iEll ++;
      }
      // evaluate at that ell
      return noise[iEll];
   }
   return 0.;
}


double int_for_C_kappaCMB(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a >1 in int_for_C_kappaCMB");

  ell       = ar[0]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res = W_kappa_cmb(a,fK)*W_kappa_cmb(a,fK)*dchi_da(a)/fK/fK;
  res = res*Pdelta(k,a);
  return res;
}

double C_kappaCMB(double l){
  double array[1] = {l};
  return int_gsl_integrate_medium_precision(int_for_C_kappaCMB,(void*)array,1./(1+redshift.shear_zdistrpar_zmin),1.0,NULL,1000);
}


/**********covariance routines ***************/

/********** var of gallens x cmblens ***************/

double cov_G_shearxkappaCMB_shearxkappaCMB(double l, double delta_l, int z1, int z3){
   double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
   C13 = C_shear_tomo_nointerp(l,z1,z3);C24 = C_kappaCMB(l);
   C14 = C_shear_CMBlensing_nointerp(l,z1);C23 = C_shear_CMBlensing_nointerp(l,z3);
   if (z1 == z3){
      N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
   }
   N24=kappa_reconstruction_noise(l);
   return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_shearxkappaCMB_shearxkappaCMB(double a,void *params)
{
  double k1,k2,fK,weights,res;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_kappa(a,fK,ar[2])*W_kappa_cmb(a,fK)*W_kappa(a,fK,ar[3])*W_kappa_cmb(a,fK)*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
  }
  res *= weights;
  return res;
}

double cov_NG_shearxkappaCMB_shearxkappaCMB(double l1,double l2, int z1, int z2){
  double a1,a2,array[4];
  int zmin;
  zmin = (z2 < z1 ? z2 : z1);
  a1 = amin_source(zmin);
  a2 = amax_source(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_shearxkappaCMB_shearxkappaCMB,(void*)array,a1,a2,NULL,1000);
}


/********** var of cmblens x projected density ***************/

double cov_clxkappaCMB_clxkappaCMB(double l, double delta_l, int n1, int n3){
   double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
   C13 = C_cl_tomo(l,n1,n3);
   C24 = C_kappaCMBxkappaCMB(l);  // !!!!!!!!!!!!!!
   C14 = C_cl_cmb(l,n1);   // !!!!!!!!!!!!!!
   C23 = C_cl_cmb(l,n3);   // !!!!!!!!!!!!!!
   if (n1 == n3){
      N13 = 1./(nlens(z1)*survey.n_gal_conversion_factor);   //!!!!!!!!!!!!!!!!!!!!!
   }
   N24=kappa_reconstruction_noise(l);
   return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_clxkappaCMB_clxkappaCMB(double a,void *params)
{
   double k1,k2,fK,weights,res;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,fK,ar[2])*W_kappa_cmb(a,fK)*W_gal(a,fK,ar[3])*W_kappa_cmb(a,fK)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_clxkappaCMB_clxkappaCMB(double l1,double l2, int n1, int n2){
   double a1,a2,array[4];
   int zmin;
   zmin = (z2 < z1 ? z2 : z1);
   a1 = amin_source(zmin);
   a2 = amax_source(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) n1;
   array[3] = (double) n2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_clxkappaCMB_clxkappaCMB,(void*)array,a1,a2,NULL,1000);
}


/********** power spectrum cmblens x cmblens ***************/
// !!!!!!!!! Could go in CMBxLSS.c
// !!!!!!!!! already coded up above

double int_for_C_kappaCMB_kappaCMB(double a, void *params)
{
   double *ar = (double *) params;
   double res,ell, fK, k;
//   if (a >= 1.0) error("a >1 in int_for_C_kappaCMB_kappaCMB");
   ell       = ar[1]+0.5;
   fK     = f_K(chi(a));
   k      = ell/fK;
   res= pow(W_kappa_cmb(a,fK), 2)*dchi_da(a)/fK/fK;
   res= res*Pdelta(k,a);
   return res;
}

double C_kappaCMB_kappaCMB(double l) {
   return int_gsl_integrate_medium_precision(int_for_C_kappaCMB_kappaCMB,(void*)array,amin_source(j),amax_source(k),NULL,1000);
}


/********** var of cmblens x cmblens ***************/

double inner_project_tri_cov_kappaCMBkappaCMB_kappaCMBkappaCMB_tomo(double a,void *params)
{
   double k1,k2,fK,weights,res;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = pow(W_kappa_cmb(a,fK), 4)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,ar[2])*pow(fK,-4.); //SSC
   }
   res *= weights;
   return res;
}

double cov_NG_kappaCMBkappaCMB_kappaCMBkappaCMB_tomo(double l1,double l2){
   double array[3];
   array[0] = l1;
   array[1] = l2;
   array[2] = survey.area/41253.0;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_kappaCMBkappaCMB_kappaCMBkappaCMB_tomo,(void*)array,NULL,1000);
}

double cov_G_kappaCMBkappaCMB_kappaCMBkappaCMB_tomo(double l){
   double C, N;
   C = C_kappaCMB_kappaCMB(l)
   N = kappa_reconstruction_noise(l)
   
   return pow((C+N), 2) / ((2.*l+1.)*delta_l*fsky);
}









