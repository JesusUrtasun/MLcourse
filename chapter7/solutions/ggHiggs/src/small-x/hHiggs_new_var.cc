#include <cmath>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>

#include <cuba.h>    // change here the CUBA path!!!!

#include "../math/integration.hh"
#include "hHiggs_new_var.hh"


using namespace std;


//const double Pi = M_PI;


dcomplex bbrn_q_ep0(dcomplex x) {
  dcomplex logx = log(x);
  dcomplex HPL2 = 0.5*logx*logx;
  return -16.*(1.+x)*(1.+x)/x*HPL2 + 32.*(x-1.)*(x-1.)/x;
}
double HiggsSmallX::factorH() {
  dcomplex ooH = 0;
  for(int i=0; i<Yi.size(); i++) {
    dcomplex Wq = 4.*Yi[i];
    dcomplex Xq = -Wq/pow(1.+sqrt(1.-Wq),2);
    ooH += bbrn_q_ep0(Xq)*Wq*Wq /64./4.;
  }
  return 1./norm(ooH);
}



// COMPLEX version using 't Hooft and Veltman's change of variable
const double PISQ6 = 1.64493406684822643647;
dcomplex CLi2(dcomplex x) {
  double x_0 = -0.30;
  double x_1 =  0.25;
  double x_2 =  0.51;
  if (x == 1.) return PISQ6;
  if (real(x) >= x_2) return PISQ6 - CLi2(1.-x) - log(x)*log(1.-x);
  if ((ABS(imag(x)) > 1.) || (real(x)*real(x) + imag(x)*imag(x) > 1.2))
    return - CLi2(1./x) - 0.5 * log(-x) * log(-x) - PISQ6 ;
  if (real(x) <= x_0) {
    dcomplex zz = log(1.-x);
    return -CLi2(-x/(1.-x)) - zz*zz/2. ;
  } else if (real(x) < x_1){
    dcomplex z = - log(1.-x);
    dcomplex temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
                   *(1.-5.*z*z/294.*(1.-7.*z*z/360.
                   *(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
                   *(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)
                    ))))))));
    return temp;
  } else return - CLi2(-x) + CLi2(x*x)/2.;
}
// GSL version of complex Li2
dcomplex gslLi2(dcomplex z) {
  gsl_sf_result re, im;
  gsl_sf_complex_dilog_e(abs(z), arg(z), &re, &im);
  return dcomplex(re.val,im.val);
}
// Li2 selector
dcomplex Li2(dcomplex z) {
  return gslLi2(z);
  //return CLi2(z);
}





dcomplex B0(double rho, dcomplex cyt) {
  double yt = real(cyt);
  dcomplex res;
  if(rho==0) res = -2.;
  else if(rho==4.*yt) res = 0.;
  else if(rho<0) {
    double sqr = sqrt(1.-4.*yt/rho);
    double sqrminusone = sqr - 1.;
    if(fabs(rho/yt)>1.e12) sqrminusone = -2.*yt/rho;
    res = -sqr * ( log(sqr+1.) - log(sqrminusone) );
  } else if(rho>4.*yt) {
    dcomplex sqr = sqrt(1.-4.*cyt/rho);
    res = -sqr * ( log(sqr+1.) - log(sqr-1.) );
  } else { // 0 < rho < 4*yt
    double sqr = sqrt((4.*yt-rho)/rho);
    res = -sqr * atan(1./sqr) *2.;
  }
  return res;
}
dcomplex dB0(double rho, dcomplex cyt) {
  double yt = real(cyt);
  dcomplex res;
  if(rho==0) res = 1./6./yt;
  else if(rho<0) {
    double sqr = sqrt(1.-4.*yt/rho);
    double sqrminusone = sqr - 1.;
    if(fabs(rho/yt)>1.e12) sqrminusone = -2.*yt/rho;
    //res = - ( ( log(sqr+1.) - log(sqr-1.) )/sqr + 1./(sqr+1.) - 1./(sqr-1.) ) *2.*yt/rho/rho;
    res = - ( log(sqr+1.) - log(sqrminusone) ) *2.*yt/rho/rho/sqr - 1./rho;
  } else if(rho>=4.*yt) {
    dcomplex sqr = sqrt(1.-4.*cyt/rho);
    //res = - ( ( log(sqr+1.) - log(sqr-1.) )/sqr + 1./(sqr+1.) - 1./(sqr-1.) ) *2.*yt/rho/rho;
    res = - ( log(sqr+1.) - log(sqr-1.) ) *2.*yt/rho/rho/sqr - 1./rho;
  } else { // 0 < rho < 4*yt
    double sqr = sqrt((4.*yt-rho)/rho);
    res = ( atan(1./sqr)/sqr -1./(sqr*sqr+1.) ) *4.*yt/rho/rho;
  }
  return res;
}
dcomplex ddB0(double rho, dcomplex cyt) {
  double yt = real(cyt);
  dcomplex res;
  if(rho==0) res = 1./30./yt/yt;
  else if(rho<0) {
    double sqr = sqrt(1.-4.*yt/rho);
    double sqrminusone = sqr - 1.;
    if(fabs(rho/yt)>1.e12) sqrminusone = -2.*yt/rho;
    res = ( log(sqr+1.) - log(sqrminusone) ) *4.*yt*(rho-3.*yt)/rho/rho/rho/(rho-4.*yt)/sqr - 2.*yt/rho/rho/rho/sqr/sqr + 1./rho/rho;
  } else if(rho>4.*yt) {
    dcomplex sqr = sqrt(1.-4.*cyt/rho);
    res = ( log(sqr+1.) - log(sqr-1.) ) *4.*yt*(rho-3.*yt)/rho/rho/rho/(rho-4.*yt)/sqr - 2.*yt/rho/rho/rho/sqr/sqr + 1./rho/rho;
  } else { // 0 < rho < 4*yt
    double sqr = sqrt((4.*yt-rho)/rho);
    res = - ( atan(1./sqr)/sqr -1./(sqr*sqr+1.) ) *8.*yt/rho/rho/rho
      - 8.*yt*yt/rho/rho/rho * ( -atan(1./sqr)/sqr/sqr/sqr -1./(sqr*sqr+1.)/sqr/sqr + 2./(sqr*sqr+1.)/(sqr*sqr+1.) );
  }
  return res;
}

double Btilde(double t, double y, dcomplex cyt) {
  //double yt = real(cyt);
  //if(t>1.e4) return -y/8./Pi/Pi;
  //if(y<0.02) return ( -1. + 4.*yt/sqrt(t*t+4.*t*yt)*atanh(sqrt(t/(t+4.*yt))) )*y/8./Pi/Pi;
  return real(B0(-t*(1.+y),cyt) - B0(-t*(1.-y),cyt));
  //dcomplex res;
  //dcomplex sqr1 = sqrt(1.+4.*yt/x1);
  //dcomplex sqr2 = sqrt(1.+4.*yt/x2);
  //res = - ( sqr1 * ( log(1.+1./sqr1) - log(1.-1./sqr1) ) - sqr2 * ( log(1.+1./sqr2) - log(1.-1./sqr2) ) ) / (16.*Pi*Pi);
  //return real(res);
}
double Btildedt(double t, double y, dcomplex cyt) {
  //double yt = real(cyt);
  return real(-(1.+y)*dB0(-t*(1.+y),cyt) + (1.-y)*dB0(-t*(1.-y),cyt));
}
double Btildedy(double t, double y, dcomplex cyt) {
  //double yt = real(cyt);
  return -t*real(dB0(-t*(1.+y),cyt) + dB0(-t*(1.-y),cyt));
}
double Btildedtdy(double t, double y, dcomplex cyt) {
  //double yt = real(cyt);
  return real(t*(ddB0(-t*(1.+y),cyt) + ddB0(-t*(1.-y),cyt))
	      - ( dB0(-t*(1.+y),cyt) +  dB0(-t*(1.-y),cyt) )
	      + y*t*( ddB0(-t*(1.+y),cyt) - ddB0(-t*(1.-y),cyt) )
	      );
}
/*
double D_0_Btilde(double x1, double x2, dcomplex cyt) {
  return -real(D_0_B0(-x1,cyt));
  double yt = real(cyt);
  dcomplex res1, res2;
  dcomplex sqr1 = sqrt(1.+4.*yt/x1);
  dcomplex sqr2 = sqrt(1.+4.*yt/x2);
  res1 = - ( ( log(1.+1./sqr1) - log(1.-1./sqr1) )/sqr1 + 1./(sqr1+1.) - 1./(sqr1-1.) ) *2.*yt/x1/x1 / (16.*Pi*Pi);
  res2 = - sqr2 * ( log(1.+1./sqr2) - log(1.-1./sqr2) ) / (16.*Pi*Pi);
  return real(res1-res2);
}
double D_1_Btilde(double x1, double x2, dcomplex cyt) {
  return real(D_0_B0(-x2,cyt));
  return -D_0_Btilde(x2,x1,cyt);
}
double D_0_1_Btilde(double x1, double x2, dcomplex cyt) {
  return 0;
  return real(D_0_B0(-x1,cyt) - D_0_B0(-x2,cyt));
  double yt = real(cyt);
  dcomplex res1, res2;
  dcomplex sqr1 = sqrt(1.+4.*yt/x1);
  dcomplex sqr2 = sqrt(1.+4.*yt/x2);
  res1 = - ( ( log(1.+1./sqr1) - log(1.-1./sqr1) )/sqr1 + 1./(sqr1+1.) - 1./(sqr1-1.) ) *2.*yt/x1/x1 / (16.*Pi*Pi);
  res2 = - ( ( log(1.+1./sqr2) - log(1.-1./sqr2) )/sqr2 + 1./(sqr2+1.) - 1./(sqr2-1.) ) *2.*yt/x2/x2 / (16.*Pi*Pi);
  return real(res1-res2);
}
*/




double Del3(double t, double y) {
  return 1.+ 4.*t + 4.*t*t*y*y;
}




//
//   y = 1
//
dcomplex C0y1(double t, dcomplex yt) {
  double sD3 = 1.+2.*t;
  dcomplex T0 = sqrt(1.-4.*yt);
  dcomplex Tp = sqrt(t +2.*yt);
  double st = sqrt(t);
  return ( - Li2(2./(1.-T0)) - Li2(2./(1.+T0)) + Li2(2.*st/(st-Tp)) + Li2(2.*st/(st+Tp)) ) /sD3;
}
dcomplex C0y1dt(double t, dcomplex yt) {
  double sD3 = 1.+2.*t;
  dcomplex Tp = sqrt(1.+2.*yt/t);
  return ( - log(1.-2./(1.+Tp))/(1.+Tp) + log(1.-2./(1.-Tp))/(1.-Tp) ) *yt/t/t/Tp/sD3 - 2.*C0y1(t,yt)/sD3;
}
dcomplex cA1y1(double t, dcomplex yt) {  // verified with Mathematica
  return C0y1(t,yt) * ( 4.*yt/(1.+2.*t) - 1. )
    + 4.*t/(1.+2.*t)/(1.+2.*t) * ( B0(-2.*t,yt) - B0(1,yt) )
    + 2./(1.+2.*t);
}
dcomplex cA1y1dt(double t, dcomplex yt) {  // verified with Mathematica
  return C0y1dt(t,yt) * ( 4.*yt/(1.+2.*t) - 1. )
    - C0y1(t,yt) * 8.*yt/(1.+2.*t)/(1.+2.*t)
    + 4.*(1.-2.*t)/(1.+2.*t)/(1.+2.*t)/(1.+2.*t) * ( B0(-2.*t,yt) - B0(1,yt) )
    - 8.*t/(1.+2.*t)/(1.+2.*t) * dB0(-2.*t,yt)
    - 4./(1.+2.*t)/(1.+2.*t);
}






/*
// Marzani's version (problematic cancellation at large t)
dcomplex kappa(double d, double xi, dcomplex yt) {
  double dp = (1.+d)/2.;
  double dm = (1.-d)/2.;
  dcomplex xp = -(xi+sqrt(xi*xi+4.*yt*xi))/2./yt;
  dcomplex xm = -(xi-sqrt(xi*xi+4.*yt*xi))/2./yt;
  return log(1-xm) * (log(1-xm*dp) - log(1-xm*dm)) + Li2(xp*dp) + Li2(xm*dp) - Li2(xp*dm) - Li2(xm*dm);
}
dcomplex C0(double t, double y, dcomplex yt) {
  double x1 = t*(1+y);
  double x2 = t*(1-y);
  double sD3 = sqrt(Del3(t,y));
  double d1 = -(1.+2.*t*y)/sD3;
  double d2 = -(1.-2.*t*y)/sD3;
  double d3 =  (1.+2.*t  )/sD3;
  return ( kappa(d1, x1, yt) + kappa(d2, x2, yt) + kappa(d3, -1., yt) ) /sD3;
}
*/
// Pasechnik-Teryaev-Szczurek's version (much more stable)
dcomplex kappa(double d, dcomplex T, double one = 1.) {
  return Li2((d-one)/(d-T)) + Li2((d-one)/(d+T)) - Li2((d+one)/(d-T)) - Li2((d+one)/(d+T));
}
dcomplex C0(double t, double y, dcomplex yt) {
  if(y==1) return C0y1(t,yt);
  double x1 = t*(1+y);
  double x2 = t*(1-y);
  double sD3 = sqrt(Del3(t,y));
  double d1 = -(1.+2.*t*y)/sD3;
  double d2 = -(1.-2.*t*y)/sD3;
  double d3 =  (1.+2.*t  )/sD3;
  dcomplex T1 = sqrt(x1+4.*yt);
  dcomplex T2 = sqrt(x2+4.*yt);
  dcomplex T3 = sqrt(1.-4.*yt);
  return ( kappa(sqrt(x1)*d1,T1,sqrt(x1)) + kappa(sqrt(x2)*d2,T2,sqrt(x2)) + kappa(d3,T3) ) /sD3;
}
//
// derivatives (considering their original form, with one=1 and T=sqrt(1+4yt/x) )
dcomplex kappadd(double d, dcomplex T, double one = 1., double *dminus1=NULL, double *dplus1=NULL, dcomplex *Tminus1=NULL) {
  double dm1 = d-one;
  if(dminus1 != NULL) dm1 = *dminus1;
  double dp1 = d+one;
  if(dplus1 != NULL) dp1 = *dplus1;
  dcomplex Tm1 = T-one;
  if(Tminus1 != NULL) Tm1 = *Tminus1;
  return one*(-one+T)/dm1/(d-T)* log((-Tm1  )/(d-T))
    +    one*(-one-T)/dm1/(d+T)* log(( T+one)/(d+T))
    -    one*( one+T)/dp1/(d-T)* log((-T-one)/(d-T))
    -    one*( one-T)/dp1/(d+T)* log(( Tm1  )/(d+T));
}
// factored out a "one"
dcomplex kappadT(double d, dcomplex T, double one = 1., dcomplex *Tminus1=NULL) {
  dcomplex Tm1 = T-one;
  if(Tminus1 != NULL) Tm1 = *Tminus1;
  return log((-Tm1  )/(d-T)) /(T-d)
    +    log(( T+one)/(d+T)) /(T+d)
    -    log((-T-one)/(d-T)) /(T-d)
    -    log(( Tm1  )/(d+T)) /(T+d);
}
dcomplex C0dt(double t, double y, dcomplex yt) {  // checked, ok
  if(y==1) return C0y1dt(t,yt);
  double x1 = t*(1+y);
  double x2 = t*(1-y);
  double sD3 = sqrt(Del3(t,y));
  double d1 = -(1.+2.*t*y)/sD3;
  double d2 = -(1.-2.*t*y)/sD3;
  double d3 =  (1.+2.*t  )/sD3;
  dcomplex T1 = sqrt(x1+4.*yt);
  dcomplex T2 = sqrt(x2+4.*yt);
  dcomplex T3 = sqrt(1.-4.*yt);
  //
  dcomplex lC0 = C0(t,y,yt);
  dcomplex res = -2.*(1.+2.*t*y*y)*lC0/sD3/sD3
    - 2.*yt/t* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 + kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /sD3;
  double d3m1 = d3-1.;
  if(t<1.e-7) d3m1 = 2.*(1.-y*y)*t*t;
  if(y>0.999) d3m1 = (1-y)*4.*t*t/(1.+2.*t)/(1.+2.*t);
  res += 2./t/sD3/sD3/sD3* kappadd(d3,T3,1,&d3m1)*2.*x1*x2 /sD3;
  if(t<1.e-5 && y>0.99) {
    double d1p1 = sqrt(x1)*( 2.*(1.-y)*t + 2.*(y*y+2.*y-3.)*t*t );
    double d2p1 = sqrt(x2)*( 2.*(1.+y)*t + 2.*(y*y-2.*y-3.)*t*t );
    res += 2./t/sD3/sD3/sD3* ( kappadd(sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)
			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y) )/sD3;
  } else
    res += 2./t/sD3/sD3/sD3* ( kappadd(sqrt(x1)*d1,T1,sqrt(x1))*x2*(1.-2.*t*y) + kappadd(sqrt(x2)*d2,T2,sqrt(x2))*x1*(1.+2.*t*y) )/sD3;
  return res;
}
dcomplex C0dy(double t, double y, dcomplex yt) {  // checked, ok
  double x1 = t*(1+y);
  double x2 = t*(1-y);
  double sD3 = sqrt(Del3(t,y));
  double d1 = -(1.+2.*t*y)/sD3;
  double d2 = -(1.-2.*t*y)/sD3;
  double d3 =  (1.+2.*t  )/sD3;
  dcomplex T1 = sqrt(x1+4.*yt);
  dcomplex T2 = sqrt(x2+4.*yt);
  dcomplex T3 = sqrt(1.-4.*yt);
  //
  dcomplex lC0 = C0(t,y,yt);
  dcomplex res = -4.*t*t*y*lC0/sD3/sD3
    - 2.*yt*t* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1/x1 - kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2/x2 ) /sD3;
  double d3m1 = d3-1.;
  if(t<1.e-7) d3m1 = 2.*(1.-y*y)*t*t;
  if(y>0.999) d3m1 = (1-y)*4.*t*t/(1.+2.*t)/(1.+2.*t);
  res += -2.*t/sD3/sD3/sD3* kappadd(d3,T3,1.,&d3m1)*2.*t*y*(1.+2.*t) /sD3;
  if(t<1.e-5 && y>0.99) {
    double d1p1 = sqrt(x1)*( 2.*(1.-y)*t + 2.*(y*y+2.*y-3.)*t*t );
    double d2p1 = sqrt(x2)*( 2.*(1.+y)*t + 2.*(y*y-2.*y-3.)*t*t );
    res += 2.*t/sD3/sD3/sD3* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*(-1.+2.*t*(y-2.))
			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*(1.+2.*t*(y+2.)) )/sD3;
  } else
    res += 2.*t/sD3/sD3/sD3* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1))*(-1.+2.*t*(y-2.))
			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2))*(1.+2.*t*(y+2.)) )/sD3;
  return res;
}
//
// second derivatives (considering their original form, with one=1 and T=sqrt(1+4yt/x) )
dcomplex kappadddd(double d, dcomplex T, double one = 1., double *dminus1=NULL, double *dplus1=NULL) {
  double dm1 = d-one;
  if(dminus1 != NULL) dm1 = *dminus1;
  double dp1 = d+one;
  if(dplus1 != NULL) dp1 = *dplus1;
  return one*one*( one-T)/dm1/pow(d-T,2)*( 1. + (2.*d-one-T)/dm1*log((-T+one)/(d-T)) )
    +    one*one*( one+T)/dm1/pow(d+T,2)*( 1. + (2.*d-one+T)/dm1*log(( T+one)/(d+T)) )
    -    one*one*(-one-T)/dp1/pow(d-T,2)*( 1. + (2.*d+one-T)/dp1*log((-T-one)/(d-T)) )
    -    one*one*(-one+T)/dp1/pow(d+T,2)*( 1. + (2.*d+one+T)/dp1*log(( T-one)/(d+T)) );
}
// factored out a "one"
dcomplex kappadTdd(double d, dcomplex T, double one = 1.) {
  return ( log((-T+one)/(d-T)) + 1. ) *one/pow(d-T,2)
    -    ( log(( T+one)/(d+T)) + 1. ) *one/pow(d+T,2)
    -    ( log((-T-one)/(d-T)) + 1. ) *one/pow(d-T,2)
    +    ( log(( T-one)/(d+T)) + 1. ) *one/pow(d+T,2);
}
// factored out a "one" squared
dcomplex kappadTdT(double d, dcomplex T, double one = 1.) {
  return ( (d-one)/( one-T) - log((-T+one)/(d-T)) ) /pow(d-T,2)
    +    ( (d-one)/( one+T) - log(( T+one)/(d+T)) ) /pow(d+T,2)
    -    ( (d+one)/(-one-T) - log((-T-one)/(d-T)) ) /pow(d-T,2)
    -    ( (d+one)/(-one+T) - log(( T-one)/(d+T)) ) /pow(d+T,2);
}
dcomplex C0dtdy(double t, double y, dcomplex yt) {
  //return gC0dtdy(t,y,yt);
  double x1 = t*(1+y);
  double x2 = t*(1-y);
  double D3 = Del3(t,y);
  double sD3 = sqrt(D3);
  double d1 = -(1.+2.*t*y)/sD3;
  double d2 = -(1.-2.*t*y)/sD3;
  double d3 =  (1.+2.*t  )/sD3;
  dcomplex T1 = sqrt(x1+4.*yt);
  dcomplex T2 = sqrt(x2+4.*yt);
  dcomplex T3 = sqrt(1.-4.*yt);
  //
  dcomplex lC0 = C0(t,y,yt);
  //dcomplex lC0dt = C0dt(t,y,yt);
  dcomplex lC0dy = C0dy(t,y,yt);
  dcomplex res = - 2.*(1.+2.*t*y*y)*lC0dy/D3 - 8.*t*y*(1.+2.*t)*lC0/D3/D3
    //
    + 8.*yt*t*y* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 + kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /D3/sD3 // deriv wrt 1/sD3
    + 2.*yt/sD3* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1/x1 - kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2/x2 ) // deriv wrt 1/x (hidden)
    - 4.*yt/D3/D3* ( kappadTdd  (sqrt(x1)*d1,T1,sqrt(x1))*(-1.+2.*t*(y-2.))/T1
    		     + kappadTdd(sqrt(x2)*d2,T2,sqrt(x2))*( 1.+2.*t*(y+2.))/T2 ) // deriv wrt d
    + 4.*yt*yt/sD3* ( (   kappadTdT(sqrt(x1)*d1,T1,sqrt(x1)) - kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 ) /T1/T1/x1
		      - ( kappadTdT(sqrt(x2)*d2,T2,sqrt(x2)) - kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /T2/T2/x2 ) // deriv wrt T
    ;
  //
  double d3m1 = d3 - 1.;
  double d1p1 = (d1 + 1.)*sqrt(x1);
  double d2p1 = (d2 + 1.)*sqrt(x2);
  if(t<1.e-7) d3m1 = 2.*(1.-y*y)*t*t;
  if(y>0.999) d3m1 = (1-y)*4.*t*t/(1.+2.*t)/(1.+2.*t);
  if(t<1.e-5 && y>0.99) {  // not mandatory
    d1p1 = sqrt(x1)*( 2.*(1.-y)*t + 2.*(y*y+2.*y-3.)*t*t );
    d2p1 = sqrt(x2)*( 2.*(1.+y)*t + 2.*(y*y-2.*y-3.)*t*t );
  }
  res += - 32.*t*y/D3/D3/D3* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)
			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y)
			       + kappadd(d3,T3,1.,&d3m1)*2.*x1*x2 ) // deriv wrt 1/D3/D3
    + 2./D3/D3* ( - kappadd(sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*(1.+2.*t*(1.-2.*y))
		  + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*(1.+2.*t*(1.+2.*y))
		  - kappadd(d3,T3,1.,&d3m1)*4.*t*y ) // deriv wrt each factors
    + 4./sD3/D3/D3/D3* ( kappadddd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)*(-1.+2.*t*(y-2.))
			 + kappadddd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y)*( 1.+2.*t*(y+2.))
			 - kappadddd(d3,T3,1.,&d3m1)*2.*x1*x2*2.*t*y*(1.+2.*t) ) // deriv wrt d
    - 4.*yt/D3/D3* ( kappadTdd  (sqrt(x1)*d1,T1,sqrt(x1))*x2*(1.-2.*t*y)/T1/x1
		     - kappadTdd(sqrt(x2)*d2,T2,sqrt(x2))*x1*(1.+2.*t*y)/T2/x2 ) // deriv wrt T
    ;
  return res;
}
// #include "hHiggs_func_new_var_P.cc"
// dcomplex C0dtdt(double t, double y, dcomplex yt) {  // done... wrong! :(
//   return gC0dtdt(t,y,yt); // use GiNaC expression for the time being
//   double x1 = t*(1+y);
//   double x2 = t*(1-y);
//   double D3 = Del3(t,y);
//   double sD3 = sqrt(D3);
//   double d1 = -(1.+2.*t*y)/sD3;
//   double d2 = -(1.-2.*t*y)/sD3;
//   double d3 =  (1.+2.*t  )/sD3;
//   dcomplex T1 = sqrt(x1+4.*yt);
//   dcomplex T2 = sqrt(x2+4.*yt);
//   dcomplex T3 = sqrt(1.-4.*yt);
//   //
//   dcomplex lC0 = C0(t,y,yt);
//   dcomplex lC0dt = C0dt(t,y,yt);
//   dcomplex res = - 2.*(1.+2.*t*y*y)*lC0dt/D3 + 8.*(2.+(4.*t-1.)*y*y+4.*t*t*y*y*y*y)*lC0/D3/D3
//     //
//     + 4.*yt/t*(1+2*t*y*y)* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 + kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /D3/sD3 // deriv wrt 1/sD3
//     + 2.*yt/t/t* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 + kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /sD3 // deriv wrt 1/t
//     + 2.*yt/t/t* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 + kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /sD3 // deriv wrt 1/x (hidden)
//     - 4.*yt/t/D3/D3* ( kappadTdd  (sqrt(x1)*d1,T1,sqrt(x1))*(y-1)*(2*t*y-1)/T1
// 		       + kappadTdd(sqrt(x2)*d2,T2,sqrt(x2))*(y+1)*(2*t*y+1)/T2 ) // deriv wrt d
//     + 4.*yt*yt/t/t/sD3* ( (   kappadTdT(sqrt(x1)*d1,T1,sqrt(x1)) - kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1 ) /T1/T1
// 			  + ( kappadTdT(sqrt(x2)*d2,T2,sqrt(x2)) - kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2 ) /T2/T2 ) // deriv wrt T
//     ;
//   //
//   double d3m1 = d3 - 1.;
//   double d1p1 = (d1 + 1.)*sqrt(x1);
//   double d2p1 = (d2 + 1.)*sqrt(x2);
//   if(t<1.e-7) d3m1 = 2.*(1.-y*y)*t*t;
//   if(y>0.999) d3m1 = (1-y)*4.*t*t/(1.+2.*t)/(1.+2.*t);
//   if(t<1.e-5 && y>0.99) {  // not mandatory
//     d1p1 = sqrt(x1)*( 2.*(1.-y)*t + 2.*(y*y+2.*y-3.)*t*t );
//     d2p1 = sqrt(x2)*( 2.*(1.+y)*t + 2.*(y*y-2.*y-3.)*t*t );
//   }
//   // res += 2./t/D3/D3* ( kappadd(sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)
//   //                    + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y)
//   //                    + kappadd(d3,T3,1,&d3m1)*2.*x1*x2 );
//   res += - 16./t/D3/D3/D3*(1+2*t*y*y)* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)
// 					 + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y)
// 					 + kappadd(d3,T3,1.,&d3m1)*2.*x1*x2 ) // deriv wrt 1/D3/D3
//     + 4./t/D3/D3* ( - kappadd(sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*y
// 		    + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*y
// 		    + kappadd(d3,T3,1.,&d3m1)*x1*x2/t ) // deriv wrt each factors
//     + 4./t/sD3/D3/D3/D3* ( kappadddd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*x2*(1.-2.*t*y)*(1-y)*(1-2*t*y)
// 			   + kappadddd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*x1*(1.+2.*t*y)*(1+y)*(1+2*t*y)
// 			   + kappadddd(d3,T3,1.,&d3m1)*2.*x1*x2 *2.*t*(1.-y*y) ) // deriv wrt d
//     - 4.*yt/t/t/D3/D3* ( kappadTdd  (sqrt(x1)*d1,T1,sqrt(x1))*x2*(1.-2.*t*y)/T1
// 		     + kappadTdd(sqrt(x2)*d2,T2,sqrt(x2))*x1*(1.+2.*t*y)/T2 ) // deriv wrt T
//     ;
//   return res;
// }
// dcomplex C0dydy(double t, double y, dcomplex yt) {  // TODO
//   return gC0dydy(t,y,yt); // use GiNaC expression for the time being
//   double x1 = t*(1+y);
//   double x2 = t*(1-y);
//   double sD3 = sqrt(Del3(t,y));
//   double d1 = -(1.+2.*t*y)/sD3;
//   double d2 = -(1.-2.*t*y)/sD3;
//   double d3 =  (1.+2.*t  )/sD3;
//   dcomplex T1 = sqrt(x1+4.*yt);
//   dcomplex T2 = sqrt(x2+4.*yt);
//   dcomplex T3 = sqrt(1.-4.*yt);
//   //
//   dcomplex lC0 = C0(t,y,yt);
//   dcomplex res = -4.*t*t*y*lC0/sD3/sD3
//     - 2.*yt*t* ( kappadT(sqrt(x1)*d1,T1,sqrt(x1))/T1/x1 - kappadT(sqrt(x2)*d2,T2,sqrt(x2))/T2/x2 ) /sD3;
//   double d3m1 = d3-1.;
//   if(t<1.e-7) d3m1 = 2.*(1.-y*y)*t*t;
//   if(y>0.999) d3m1 = (1-y)*4.*t*t/(1.+2.*t)/(1.+2.*t);
//   res += -2.*t/sD3/sD3/sD3* kappadd(d3,T3,1.,&d3m1)*2.*t*y*(1.+2.*t) /sD3;
//   if(t<1.e-5 && y>0.99) {
//     double d1p1 = sqrt(x1)*( 2.*(1.-y)*t + 2.*(y*y+2.*y-3.)*t*t );
//     double d2p1 = sqrt(x2)*( 2.*(1.+y)*t + 2.*(y*y-2.*y-3.)*t*t );
//     res += 2.*t/sD3/sD3/sD3* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1),NULL,&d1p1)*(-1.+2.*t*(y-2.))
// 			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2),NULL,&d2p1)*(1.+2.*t*(y+2.)) )/sD3;
//   } else
//     res += 2.*t/sD3/sD3/sD3* ( kappadd  (sqrt(x1)*d1,T1,sqrt(x1))*(-1.+2.*t*(y-2.))
// 			       + kappadd(sqrt(x2)*d2,T2,sqrt(x2))*(1.+2.*t*(y+2.)) )/sD3;
//   return res;
// }
//
//
//
//
dcomplex cA1(double t, double y, dcomplex yt) {  // checked in y = 0, 1
  double D3 = Del3(t,y);
  return C0(t,y,yt)/D3 * ( 4.*yt*(1.+2.*t) - (1.+2.*t)*(1.+2.*t) + 12.*(1.+2.*t)*t*t*(1.-y*y)/D3 )
    + 2./D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( t - 6.*t*t*(1.-y*y)/D3 )
    + 2.*t*y/D3 * Btilde(t,y,yt) * ( 1. + 12.*t*t*(1.-y*y)/D3 )
    + 2.*(1.+2.*t)/D3;
}
dcomplex cA1dt(double t, double y, dcomplex yt) {  // checked in y = 0, 1 (even if in y=1 C0dt fails)
  double D3 = Del3(t,y);
  dcomplex lC0 = C0(t,y,yt);
  dcomplex lC0dt = C0dt(t,y,yt);
  return ( lC0dt/D3 - 4.*lC0*(1.+2.*t*y*y)/D3/D3 )
    * ( 4.*yt*(1.+2.*t) - (1.+2.*t)*(1.+2.*t) + 12.*(1.+2.*t)*t*t*(1.-y*y)/D3 )
    + lC0/D3 * ( 8.*yt - 4.*(1.+2.*t) + 24.*t*(1.+5.*t+8.*t*t+4.*t*t*t*y*y)*(1.-y*y)/D3/D3 )
    - 8.*(1.+2.*t*y*y)/D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( t - 6.*t*t*(1.-y*y)/D3 )
    + 2./D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( 1. - 12.*t*(1.+2.*t)*(1.-y*y)/D3/D3 )
    - 2./D3 * ( dB0(-t*(1.+y),yt) + dB0(-t*(1.-y),yt) + y*(dB0(-t*(1.+y),yt)-dB0(-t*(1.-y),yt)) ) * ( t - 6.*t*t*(1.-y*y)/D3 )
    + 2.*y*(1.-4.*t*t*y*y)/D3/D3 * Btilde(t,y,yt) * ( 1 + 12.*t*t*(1.-y*y)/D3 )
    + 2.*t*y/D3 * Btilde(t,y,yt) * 24.*t*(1.+2.*t)*(1.-y*y)/D3/D3
    + 2.*t*y/D3 * Btildedt(t,y,yt) * ( 1 + 12.*t*t*(1.-y*y)/D3 )
    - 4.*(1.+4.*t*y*y*(1.+t))/D3/D3;
}
dcomplex cA1dy(double t, double y, dcomplex yt) {  // checked in y = 0, 0.5, 1 (even if in y=1 C0dy fails)
  if(y<0.01 && t>1.e4) return 0.;  ///////
  double D3 = Del3(t,y);
  dcomplex lC0 = C0(t,y,yt);
  return ( C0dy(t,y,yt)/D3/D3 - lC0*16.*t*t*y/D3/D3/D3 )
    * ( 4.*yt*(1.+2.*t)*D3 - (1.+2.*t)*(1.+2.*t)*D3 + 12.*(1.+2.*t)*t*t*(1.-y*y) )
    + lC0/D3/D3 * ( 8.*t*t*y*(4.*yt*(1.+2.*t) - (1.+2.*t)*(1.+2.*t)) - 24.*(1.+2.*t)*t*t*y )
    - 32.*t*t*y/D3/D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( t*D3 - 6.*t*t*(1.-y*y) )
    - 2.*t/D3 * ( dB0(-t*(1.+y),yt) - dB0(-t*(1.-y),yt) ) * ( t - 6.*t*t*(1.-y*y)/D3 )
    + 2./D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( 8.*t*t*t*y + 12.*t*t*y )
    + 2.*t*(1.+4.*t-12.*t*t*y*y)/D3/D3/D3 * Btilde(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*t*y/D3/D3 * Btildedy(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*t*y/D3/D3 * Btilde(t,y,yt) * ( 8.*t*t*y - 24.*t*t*y )
    - 16.*t*t*y*(1.+2.*t)/D3/D3
    ;
}
dcomplex cA1dtdy(double t, double y, dcomplex yt) {  // cross-checked here and with Mathematica 
  if(y<0.01 && t>1.e4) return 0.;  ///////
  double D3 = Del3(t,y);
  dcomplex lC0     = C0  (t,y,yt);
  dcomplex lC0dt   = C0dt(t,y,yt);
  dcomplex lC0dy   = C0dy(t,y,yt);
  dcomplex lC0dtdy = C0dtdy(t,y,yt);
  return ( lC0dtdy/D3/D3 - lC0dt*16.*t*t*y/D3/D3/D3
	   - 4.*lC0dy*(1.+2.*t*y*y)/D3/D3/D3 + 16.*t*y*lC0*(-1.+2.*t+8.*t*t*y*y)/D3/D3/D3/D3 )
    * ( 4.*yt*(1.+2.*t)*D3 - (1.+2.*t)*(1.+2.*t)*D3 + 12.*(1.+2.*t)*t*t*(1.-y*y) )
    + ( lC0dt/D3/D3 - 4.*lC0*(1.+2.*t*y*y)/D3/D3/D3 )
    * ( 8.*t*t*y*(4.*yt*(1.+2.*t) - (1.+2.*t)*(1.+2.*t)) - 24.*(1.+2.*t)*t*t*y )
    //
    + ( lC0dy/D3/D3/D3 - lC0*24.*t*t*y/D3/D3/D3/D3 )
    * ( (8.*yt - 4.*(1.+2.*t))*D3*D3 + 24.*t*(1.+5.*t+8.*t*t+4.*t*t*t*y*y)*(1.-y*y) )
    + lC0/D3/D3/D3 * ( (8.*yt - 4.*(1.+2.*t))*16.*t*t*y*D3 + 24.*t*8.*t*t*t*y*(1.-y*y) - 48.*t*(1.+5.*t+8.*t*t+4.*t*t*t*y*y)*y )
    //
    - 8.*( 4.*t*y/D3/D3/D3 - 24.*t*t*y*(1.+2.*t*y*y)/D3/D3/D3/D3 )
    * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( t*D3 - 6.*t*t*(1.-y*y) )
    + 8.*(1.+2.*t*y*y)/D3/D3/D3 * t*( dB0(-t*(1.+y),yt) - dB0(-t*(1.-y),yt) ) * ( t*D3 - 6.*t*t*(1.-y*y) )
    - 8.*(1.+2.*t*y*y)/D3/D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( 8.*t*t*t*y + 12.*t*t*y )
    //
    - 48.*t*t*y/D3/D3/D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( D3*D3 - 12.*t*(1.+2.*t)*(1.-y*y) )
    - 2./D3/D3/D3 * t*( dB0(-t*(1.+y),yt) - dB0(-t*(1.-y),yt) ) * ( D3*D3 - 12.*t*(1.+2.*t)*(1.-y*y) )
    + 2./D3/D3/D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( 16.*t*t*y*D3 + 24.*t*(1.+2.*t)*y )
    //
    + 32.*t*t*y/D3/D3/D3 * ( dB0(-t*(1.+y),yt) + dB0(-t*(1.-y),yt) + y*(dB0(-t*(1.+y),yt)-dB0(-t*(1.-y),yt)) )
    * ( t*D3 - 6.*t*t*(1.-y*y) )
    - 2./D3/D3 * ( -t*(ddB0(-t*(1.+y),yt) - ddB0(-t*(1.-y),yt))
		   + dB0(-t*(1.+y),yt) - dB0(-t*(1.-y),yt)
		   - y*t*( ddB0(-t*(1.+y),yt) + ddB0(-t*(1.-y),yt) )
		   ) * ( t*D3 - 6.*t*t*(1.-y*y) )
    - 2./D3/D3 * ( dB0(-t*(1.+y),yt) + dB0(-t*(1.-y),yt) + y*(dB0(-t*(1.+y),yt)-dB0(-t*(1.-y),yt)) ) * ( 8.*t*t*t*y + 12.*t*t*y )
    //
    + ( 2.*(1.-4.*t*t*y*y)/D3/D3/D3 - 16.*t*t*y*y/D3/D3/D3 - 48.*t*t*y*y*(1.-4.*t*t*y*y)/D3/D3/D3/D3 )
    * Btilde(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*y*(1.-4.*t*t*y*y)/D3/D3/D3 * Btildedy(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*y*(1.-4.*t*t*y*y)/D3/D3/D3 * Btilde(t,y,yt) * ( 8.*t*t*y - 24.*t*t*y )
    //
    + 2.*t/D3/D3/D3 * Btilde(t,y,yt) * 24.*t*(1.+2.*t)*(1.-y*y)
    - 48.*t*t*t*y*y/D3/D3/D3/D3 * Btilde(t,y,yt) * 24.*t*(1.+2.*t)*(1.-y*y)
    + 2.*t*y/D3/D3/D3 * Btildedy(t,y,yt) * 24.*t*(1.+2.*t)*(1.-y*y)
    - 2.*t*y/D3/D3/D3 * Btilde(t,y,yt) * 48.*t*(1.+2.*t)*y
    //
    + 2.*t/D3/D3 * Btildedt(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    - 32.*t*t*t*y*y/D3/D3/D3 * Btildedt(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*t*y/D3/D3 * Btildedtdy(t,y,yt) * ( D3 + 12.*t*t*(1.-y*y) )
    + 2.*t*y/D3/D3 * Btildedt(t,y,yt) * ( 8.*t*t*y - 24.*t*t*y )
    //
    - 32.*t*y*(1.+t)/D3/D3
    + 64.*t*t*y*(1.+4.*t*y*y*(1.+t))/D3/D3/D3;
}
//
//
dcomplex cA3(double t, double y, dcomplex yt) {  // verified with Mathematica
  double D3 = Del3(t,y);
  return C0(t,y,yt)/D3 * ( 8.*yt - 4. - 4.*t + 6.*(1.+2.*t)*(1.+2.*t)/D3 )
    + 2./D3 * ( B0(-t*(1.+y),yt) + B0(-t*(1.-y),yt) - 2.*B0(1,yt) ) * ( 1. - 3.*(1.+2.*t)/D3 )
    + 12.*t*y * Btilde(t,y,yt) * (1.+2.*t)/D3/D3
    + 4./D3;
}








//
//   y = 0
//
dcomplex C0y0(double t, dcomplex yt) {
  return C0(t,0,yt);
}
dcomplex C0y0dt(double t, dcomplex yt) {
  //return C0dt(t,0,yt);
  double sD3 = sqrt(Del3(t,0));
  double d1 = -1./sD3;
  double d3 = (1.+2.*t)/sD3;
  dcomplex T1 = sqrt(t +4.*yt);
  dcomplex T3 = sqrt(1.-4.*yt);
  //
  dcomplex lC0y0 = C0y0(t,yt);
  dcomplex res = -2.*lC0y0/sD3/sD3;
  if(t<1.e-7) {
    double d3m1 = 2.*t*t;
    res += 4./sD3/sD3/sD3* ( kappadd(sqrt(t)*d1,T1,sqrt(t)) + kappadd(d3,T3,1,&d3m1)*t ) /sD3
      - 4.*yt/t* kappadT(sqrt(t)*d1,T1,sqrt(t))/T1 /sD3;
  } else if(t > 1.e12) { // not used yet..
    dcomplex T1m1 = 2.*yt/t;
    res += 4./sD3/sD3/sD3* ( kappadd(d1,T1,1,NULL,NULL,&T1m1) + kappadd(d3,T3)*t ) /sD3
      - 4.*yt/t* kappadT(d1,T1,1,&T1m1)/sqrt(t)/T1 /sD3;
  } else
    res += 4./sD3/sD3/sD3* ( kappadd(sqrt(t)*d1,T1,sqrt(t)) + kappadd(d3,T3)*t ) /sD3
      - 4.*yt/t* kappadT(sqrt(t)*d1,T1,sqrt(t))/T1 /sD3;
  return res;
}
dcomplex cA1y0(double t, dcomplex yt) {  // verified with Mathematica
  double D3 = 1+4.*t;
  return C0y0(t,yt)/D3 * ( 4.*yt*(1.+2.*t) - 1. - 4.*t - 4.*t*t*(1.- 3.*(1.+2.*t)/D3 ) )
    + 4./D3 * ( B0(-t,yt) - B0(1,yt) ) * ( t - 6.*t*t/D3 )
    + 2.*(1.+2.*t)/D3;
}
dcomplex cA1y0dt(double t, dcomplex yt) {  // verified with Mathematica
  double D3 = 1+4.*t;
  dcomplex lC0y0   = C0y0  (t,yt);
  dcomplex lC0y0dt = C0y0dt(t,yt);
  return ( lC0y0dt/D3 - 4.*lC0y0/D3/D3 ) * ( 4.*yt*(1.+2.*t) - 1. - 4.*t - 4.*t*t*(1.- 3.*(1.+2.*t)/D3 ) )
    + lC0y0/D3 * ( 8.*yt - 4. - 8.*t*(1.- 3.*(1.+2.*t)/D3 ) - 24.*t*t/D3/D3 )
    + 4./D3 * ( -dB0(-t,yt) - 4.*(B0(-t,yt) - B0(1,yt))/D3 ) * ( t - 6.*t*t/D3 )
    + 4./D3 * ( B0(-t,yt) - B0(1,yt) ) * ( 1. - 12.*t/D3 + 24.*t*t/D3/D3 )
    - 4./D3/D3;
}








// general case for h_k
double f1y0dt(double t, vector<dcomplex> Y) {
  dcomplex A1y0   = 0;
  dcomplex A1y0dt = 0;
  for(int i=0; i<Y.size(); i++) {
    A1y0   += Y[i]*cA1y0  (t,Y[i]);
    A1y0dt += Y[i]*cA1y0dt(t,Y[i]);
  }
  return 2.*real(A1y0*conj(A1y0dt));
}
double f1dt(double t, double y, vector<dcomplex> Y) {
  dcomplex A1   = 0;
  dcomplex A1dt = 0;
  for(int i=0; i<Y.size(); i++) {
    A1   += Y[i]*cA1  (t,y,Y[i]);
    A1dt += Y[i]*cA1dt(t,y,Y[i]);
  }
  return 2.*real(A1*conj(A1dt));
}
double f1dtdy(double t, double y, vector<dcomplex> Y) {
  dcomplex A1     = 0;
  dcomplex A1dt   = 0;
  dcomplex A1dy   = 0;
  dcomplex A1dtdy = 0;
  for(int i=0; i<Y.size(); i++) {
    A1     += Y[i]*cA1    (t,y,Y[i]);
    A1dt   += Y[i]*cA1dt  (t,y,Y[i]);
    A1dy   += Y[i]*cA1dy  (t,y,Y[i]);
    A1dtdy += Y[i]*cA1dtdy(t,y,Y[i]);
  }
  return 2.*real(A1dt*conj(A1dy)) + 2.*real(A1*conj(A1dtdy));
}
double f2(double t, double y, vector<dcomplex> Y) {
  dcomplex dA3 = 0;
  for(int i=0; i<Y.size(); i++) {
    dA3 += Y[i]*cA3(t,y,Y[i]);
  }
  return 2.*norm(dA3); // = 2*|A3|^2
}

//
// 2-dimensional integrand
// new simplified verision (must use only2d=true)
double integrand2d(double t, double y, int kp, int km, vector<dcomplex> Y, double lFQ) {
  double lf1dtdy  = f1dtdy(t,y,Y);
  double lf1dt    = f1dt  (t,y,Y);
  double lyf1dtdy = lf1dt + y*lf1dtdy;
  double res;
  double lf2 = f2(t,y,Y);
  double Lp = log(t*t*(1.-y*y)) - 2*lFQ;
  double Lm = log((1.+y)/(1.-y));
  res = pow(Lp,kp) * pow(Lm,km) * ( - lyf1dtdy
				    - 2*kp /Lp * lf1dt
				    + 4.*kp*(kp-1.)*t /Lp/Lp * lf2
				    - 4.*km*(km-1.)*t /Lm/Lm * lf2
				    - km/(1.+kp) * Lp/Lm * lf1dtdy
				    );
  return res/gsl_sf_fact(kp)/gsl_sf_fact(km);
}
// 1-dimensional integral
// double integrand1d(double t, int kp) {
//   if(t>1.e10) return 0; ///////////
//   double lf1y0dt = f1y0dt(t);
//   double res = -2.*prefH * pow(2.*log(t),kp) * lf1y0dt / gsl_sf_fact(kp);
//   if(res!=res) return 0; ////////// if res is NaN, then return 0
//   return res;
// }


// parameters for generalized GSL change of variable
const double pbeta1 = 2.4;
const double pbeta2 = 3.2;
//
// 1-dimensional integrand
// double integrand1d_ch_var(double x, void *p) {
//   int kp = *(int*)p;
//   // generalized GSL change of variables
//   double t = pow(1.-x,pbeta2)/pow(x,pbeta1);
//   double jac = (pbeta1+x*(pbeta2-pbeta1)) * pow(1.-x,pbeta2-1.) / pow(x,pbeta1+1.);
//   double res = integrand1d(t,kp);
//   return res*jac;
// }
double integrand2d_ch_var(double x, double y, int kp, int km, vector<dcomplex> Y, double lFQ) {
  // generalized GSL change of variables
  double t = pow(1.-x,pbeta2)/pow(x,pbeta1);
  double jac = (pbeta1+x*(pbeta2-pbeta1)) * pow(1.-x,pbeta2-1.) / pow(x,pbeta1+1.);
  //
  //if(t>1.e10) return 0; ///////////
  //
  // y change of variable: yold = 1-y^2
  double yold = 1.-y*y;
  jac *= 2.*y;
  //
  double res = integrand2d(t,yold,kp,km,Y,lFQ)*jac;
  //
  //if(res!=res) cout << t << "  " << yold << "  " << res << endl;
  if(res!=res) return 0; ////////// if res is NaN, then return 0
  //cout << x << "   \t" << t << "   \t" << res << "   \t" << jac << endl;
  //
  return res;
}
struct kpars {int kp, km; vector<dcomplex> Y; double lFQ; };
int integrandCUBA(const int *ndim, const double x[], const int *ncomp, double res[], void *pars) {
  kpars k = *(kpars*)pars;
  res[0] = integrand2d_ch_var(x[0],x[1],k.kp,k.km,k.Y,k.lFQ);
  return 0;
}

// return the coefficient of the term Mp^kp Mm^km
double HiggsSmallX::h_kp_km(int kp, int km, double *err) {
  double res = 0;
  kpars k = {kp, km, Yi, 2*log(_muFrat)};
  //
  double cres[1], error[1], prob[1];
  double epsrel=5.e-4, epsabs=1.e-10;
  int last = 4;
  int verbose = 0;
  int nregions, neval, fail;
  //system("export CUBACORES=8");
  Cuhre(2, 1, integrandCUBA, &k, 1,
   	epsrel, epsabs, verbose | last,
   	0, 500000, 9,
	NULL, NULL,
  	&nregions, &neval, &fail, cres, error, prob);
  // Suave(2, 1, integrandCUBA, &k,
  // 	epsrel, epsabs, verbose | last, 0,
  // 	0, 500000, 1000, 25.,
  // 	&nregions, &neval, &fail, cres, error, prob);
  // Vegas(2, 1, integrandCUBA, &k,
  //   	epsrel, epsabs, verbose | last, 0,
  //   	0, 500000, 1000, 500, 1000,
  //   	0, NULL,
  //   	&neval, &fail, cres, error, prob);
  //  
  if(fail == -1) {
    cout << "  ERROR: dimension out of range." << endl;
  } else if(fail>0) {
    cout << "  ERROR: precision not reached." << endl;
  }
  res += cres[0];
  if(err != NULL) *err = prefH*error[0];
  return prefH*res;
}







// special case for h_1
double f1y1dt(double t, vector<dcomplex> Y) {
  dcomplex A1y1 = 0;
  dcomplex A1y1dt = 0;
  for(int i=0; i<Y.size(); i++) {
    A1y1   += Y[i]*cA1y1  (t,Y[i]);
    A1y1dt += Y[i]*cA1y1dt(t,Y[i]);
  }
  return 2.*real(A1y1*conj(A1y1dt));
}
double integrand_h1(double x, void *par) {
  vector<dcomplex> Y = *(vector<dcomplex> *)par;
  double t = (1.-x)/x;
  double df10 = f1y1dt(t,Y);
  double res = -2.*log(2.*t)*df10;
  return res/x/x;
}
double HiggsSmallX::h1() {
  double prec = 1e-6, error;
  double res = integration::gauss(integrand_h1, 0., 1., prec, &error, &Yi);
  return prefH*res - 4*log(_muFrat);
}





// initialization
void HiggsSmallX::init(double mH, double mt, double mb, double mc, double muFrat) {
  _muFrat = muFrat;
  //
  double yt = pow(mt/mH,2.);
  double yb = pow(mb/mH,2.);
  double yc = pow(mc/mH,2.);
  //
  dcomplex YT = yt;
  if(yt < 0.25) YT += I*1.e-12;  // is it really needed???
  dcomplex YB = yb;
  if(yb < 0.25) YB += I*1.e-12;  // is it really needed???
  dcomplex YC = yc;
  if(yc < 0.25) YC += I*1.e-12;  // is it really needed???
  //
  Yi.clear();
  if(mt>0) Yi.push_back(YT);
  if(mb>0) Yi.push_back(YB);
  if(mc>0) Yi.push_back(YC);
  //
  prefH = factorH();
}

double conversion(int k1, int k2, int km) {
  if(km % 2 != 0) return 0;
  double k12 = k1+k2;
  double res = 0;
  for(int b=0; b<=km; b++) {
    if(b<km-k2) continue;
    if(b>k1)    continue;
    res += pow(-1.,b)*gsl_sf_choose(km,b)*gsl_sf_choose(k12-km,k1-b);
  }
  return res/pow(2.,k12);
}

double h_k1_k2(int k1, int k2, double **hkpm) {
  double res = 0;
  for(int km=0; km<=k1+k2; km+=2) {
    //cout << conversion(k1, k2, km) << "  " << hkpm[k1+k2-km][km]<< endl;
    //cout << km << "  " << k1 << "  " << k2 << endl;
    //cout << hkpm[k1+k2-km][km]<< endl;
    res += hkpm[k1+k2-km][km] * conversion(k1, k2, km);
  }
  return res;
}

vector<vector<double> > HiggsSmallX::hCoeffs(int maxorder) {
  int Kmax = maxorder+1;
  double **hkpm = new double* [Kmax];
  for(int i=0; i<Kmax; i++) hkpm[i] = new double[Kmax];
  double error;
  for(int kp=0; kp<Kmax; kp++) {
    for(int km=0; km<Kmax-kp; km+=2) {
      hkpm[kp][km] = h_kp_km(kp, km, &error);
      //cout << kp << "   " << km << "   " << hkpm[kp][km] << "  ( " << error/hkpm[kp][km] << " relative error )" << endl;
    }
  }
  vector<vector<double> > hk;
  hk.resize(Kmax);
  for(int k12=0; k12<Kmax; k12++) {
    hk[k12].resize(k12+1);
    for(int k1=0; k1<=k12; k1++) {
      int k2 = k12-k1;
      //if(k2>k1) continue;
      //hk12[k1][k2] = h_k1_k2(k1, k2, hkpm);
      hk[k12][k2] = h_k1_k2(k1, k2, hkpm);
      //cout << k1 << "   " << k2 << "   " << hk12[k1][k2] << endl;
    }
  }
  return hk;
}
