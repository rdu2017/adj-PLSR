

README for adj PLSR method
=====================

This file serves as an breif instruction for using provided R codes to carry out the adjusted PLSR approach to utilize additional
exposure information in study the assoication to an health outcome, and to simulate a dataset using the assumed bilinear modle 
equations in a PLSR. 

cite: An adjusted partial least squares regression framework to utilize additional exposure information in environmental mixture data analysis
contact: rdu@uams.edu


Description of the R codes: 
=============================================

1. _adjPLS1 appr.R: contains two R functions.
   adjPLS1: carry out the model fitting of the adj PLSR approach for data with univariate outcome
   stat.adjPLS1.wiBT: compute the test statistic of the adj PLSR approach for data with univariate outcome

2. _adjPLS2 appr.R: contains two R functions.
   adjPLS2: carry out the model fitting of the adj PLSR approach for data with multivariate outcome
   stat.adjPLS2.wiBT: compute the test statistic of the adj PLSR approach for data with multivariate outcome

3. _simu func.R: contains one R function.
   makdata: simulate a dataset using the bilinear modle equations assumed in PLSR
  (note: this R code requires to load and attach the 'MASS' package into the environment, so one needs to install the package is it's not.)

4. _oth funcs.R: ontains other useful R functions. 
   round2: round up a value to specified place 
   rmse: calculate the root mean squared error
   pval.2side: calculate a 2-sided p-value following the resampling test procedure
   PLS1: an auxiliary function to fit a gneral PLSR

5. _example.R: provides an example for using above R functions to implement the adj. PLSR approach on simulated datasets. 



Created June 26, 2020 
