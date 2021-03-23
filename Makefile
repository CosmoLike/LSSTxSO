cfastcov_dir := cfastcov/
cfastcov := $(cfastcov_dir)twobessel.c $(cfastcov_dir)utils.c $(cfastcov_dir)utils_complex.c
opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 \
-std=gnu99 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass
opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
opt_puma := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.6/include -L/cm/shared/uaapps/gsl/2.6/lib \
-lfftw3 -lgsl -lgslcblas -lm -g -O3 \
-ffast-math -funroll-loops -std=gnu99 -L../cosmolike_core/class -lclass
cfftlog_dir := ../cosmolike_core/cfftlog/
cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c


home: 
	make home_shared
	make home_datav
	make home_cov

home_datav:
	gcc like_fourier.c -o ./like_fourier $(opt_home)

home_datav_wfirstxso:
	gcc like_fourier_wfirstxso.c -o ./like_fourier_wfirstxso $(opt_home) -DNOMPP


home_datav_1sample:
	gcc like_fourier.c -o ./like_fourier_1sample $(opt_home) -DONESAMPLE

home_datav_1sample_outlier:
	gcc like_fourier_fast_outlier.c -o ./like_fourier_1sample_outlier $(opt_home) -DONESAMPLE
home_datav_1sample_nooutlier:
	gcc like_fourier_fast_nooutlier.c -o ./like_fourier_1sample_nooutlier $(opt_home) -DONESAMPLE

home_des:
	gcc like_fourier_desxplanck.c -o ./like_fourier_desxplanck $(opt_home)
	
home_deslib:
	gcc -shared -o like_fourier_desxplanck.so -fPIC like_fourier_desxplanck.c $(opt_home)

home_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_home)

home_cov_wfirstxso:
	gcc compute_covariances_fourier_wfirstxso.c -o ./compute_covariances_fourier_wfirstxso $(opt_home)

home_cov_binned:
	gcc compute_covariances_fourier_binned.c -o ./compute_covariances_fourier_binned $(opt_home)


home_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_home) -DONESAMPLE
home_cov_1sample_ztrue:
	gcc compute_covariances_fourier_ztrue.c -o ./compute_covariances_fourier_1sample_ztrue $(opt_home) -DONESAMPLE


home_shared:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(opt_home)
	gcc -shared -o like_fourier_1sample.so -fPIC like_fourier.c $(opt_home) -DONESAMPLE

home_shared_wfirstxso:
	gcc -shared -o like_fourier_wfirstxso.so -fPIC like_fourier_wfirstxso.c $(opt_home)

home_shared_fast:
	gcc -shared -o like_fourier.so -fPIC like_fourier_fast.c $(opt_home)
	gcc -shared -o like_fourier_1sample.so -fPIC like_fourier_fast.c $(opt_home) -DONESAMPLE
home_datav_fast:
	gcc like_fourier_fast.c -o ./like_fourier $(opt_home)


ocelote:
	make ocelote_shared 
	make ocelote_datav
	make ocelote_cov

ocelote_datav:
	gcc like_fourier.c -o ./like_fourier $(opt_ocelote)

ocelote_cov:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier $(opt_ocelote)
ocelote_cov_binned:
	gcc compute_covariances_fourier_binned.c -o ./compute_covariances_fourier_binned $(opt_ocelote)

ocelote_cov_wfirstxso:
	gcc compute_covariances_fourier_wfirstxso.c -o ./compute_covariances_fourier_wfirstxso $(opt_ocelote)

ocelote_cov_1sample:
	gcc compute_covariances_fourier.c -o ./compute_covariances_fourier_1sample $(opt_ocelote) -DONESAMPLE


ocelote_shared:
	gcc -shared -o like_fourier.so -fPIC like_fourier.c $(opt_ocelote)
	gcc -shared -o like_fourier_1sample.so -fPIC like_fourier.c $(opt_ocelote) -DONESAMPLE

ocelote_shared_fast:
	gcc -shared -o like_fourier.so -fPIC like_fourier_fast.c $(opt_ocelote)
	gcc -shared -o like_fourier_1sample.so -fPIC like_fourier_fast.c $(opt_ocelote) -DONESAMPLE

ocelote_shared_wfirstxso:
	gcc -shared -o like_fourier_wfirstxso.so -fPIC like_fourier_wfirstxso.c $(opt_ocelote)


ocelote_des:
	gcc like_fourier_desxplanck.c -o ./like_fourier_desxplanck $(opt_ocelote)
	
ocelote_deslib:
	gcc -shared -o like_fourier_desxplanck.so -fPIC like_fourier_desxplanck.c $(opt_ocelote)


###### Puma

puma_shared_fast:
	gcc -shared -o like_fourier.so -fPIC like_fourier_fast.c $(opt_puma)
	gcc -shared -o like_fourier_1sample.so -fPIC like_fourier_fast.c $(opt_puma) -DONESAMPLE
