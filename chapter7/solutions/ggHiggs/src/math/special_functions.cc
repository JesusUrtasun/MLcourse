#include "special_functions.hh"
#include <cmath>
#include <iostream>
#include <gsl/gsl_sf.h>
#include "../parameters.hh"


double Li2(double x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1.) return ZETA2;
 if (x <= x_0){ 
   double temp = log(fabs(1.0-x));
   return -Li2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double z = - log(1.0-x);
   double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
   return temp; }
   else if (x < x_2) return - Li2(-x) + Li2(x*x)/2.0 ;
   else { return ZETA2 - Li2(1.0-x) 
                  - log(fabs(x))*log(fabs(1.0-x)) ; }
}
double Li3(double x){
 double x_0 = -1.0;
 double x_1 = -0.85;
 double x_2 = 0.25;
 double x_3 = 0.63;
 double x_4 =  1.0;
 if (x == 1.) return ZETA3;
 if (x == -1.) return - 0.75 * ZETA3;
 if (x <= x_0){ 
   double lnx = log(-x);
   return Li3(1.0/x) - ZETA2*lnx - lnx*lnx*lnx/6.0; }
 else if (x < x_1){
   return Li3(x*x)/4.0 - Li3(-x); }
   else if (x < x_2){
     double z = - log(1.0-x);
     double temp = z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
                    *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
                    *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
                    *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
                    *(1.0-53598548.0*z/524808375.0
                    *(1.0+22232925.0*z/107197096.0
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return Li3(x*x)/4.0 - Li3(-x); }
       else if (x < x_4){
         double ln1x = log(1.0-x); 
         return -Li3(1.0-x) - Li3(-x/(1.0-x)) + ZETA3 + ZETA2*ln1x
	   - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0; }
       else { 
         double lnx = log(x);
         return Li3(1./x) + 2.0*ZETA2*lnx - lnx*lnx*lnx/6.0; }
}


// generalisation of binomial to real numbers
double binom(double n, double m) {
  if(n<0 && -n==int(-n)) return pow(-1.,m) * binom(-n+m-1,m);
  return gsl_sf_gamma(n+1)/gsl_sf_gamma(m+1)/gsl_sf_gamma(n-m+1);
}

void reCoeff(double *c_in, double *c_out, double power, int truncation) {
  if(power == 0) {
    for(int i=0; i<=truncation; i++) c_out[i] = c_in[i];
    return;
  }
  for(int i=0; i<=truncation; i++) {
    c_out[i] = 0;
    for(int j=0; j<=i; j++) {
      if(power>0 && power==int(power) && j>power) {
	continue;
      }
      //cout << setw(3) << i << setw(3) << j << setw(12) << c_in[i-j] << setw(12) << binom(power,j) << endl;
      c_out[i] += c_in[i-j] * pow(-1.,j) * binom(power,j);
    }
    //cout << c_out[i] << endl;
  }
  //cout << endl;
}

