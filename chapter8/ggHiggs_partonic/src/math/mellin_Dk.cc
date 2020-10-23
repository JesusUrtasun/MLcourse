/*--------------------------

The aim of this library is the computation of the Mellin transforms of

Dk(z)  = [ log^k (1-z) / (1-z) ]_+

DMk(z) = [ log^k (log(1/z)) / log(1/z) ]_+

DBk(z) = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+'

where the prime (') in the last denotes that, in fact,
the log sqrt(z) term is outside the plus-distribution.

---------------------------*/


#include <iostream>
#include <cmath>
#include <gsl/gsl_sf.h>
#include "mellin_Dk.hh"





dcomplex psi0(dcomplex Z){
  dcomplex SUB = 0. ;
  dcomplex ZZ = Z;
  if(ABS(imag(ZZ))<10.) { // if too close to the real axis...
  label1:
    if(real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
      SUB = SUB - 1./ ZZ;
      ZZ = ZZ + 1.;
      goto label1;
    }
  }
  dcomplex RZ = 1./ ZZ;
  dcomplex DZ = RZ * RZ;
  // SUB + asympt expansion (Abramowitz, Stengun, 6.3.18)
  return SUB + log(ZZ) - 0.5 * RZ - DZ/5040. * ( 420.+ DZ * ( - 42. + DZ * (20. - 21. * DZ) ) );
}

dcomplex psi_n(int M, dcomplex Z){
  if(M==0) return psi0(Z);
  int K1, K2;
  dcomplex SUB = 0. , SUBM;
  dcomplex ZZ = Z;
  if(ABS(imag(ZZ))<10.) { // if too close to the real axis...
  label2:
    SUBM = -1./ZZ;
    for(K1=1; K1<=M; K1++) {
      SUBM = - SUBM * K1 / ZZ;
    }
    if(real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
      SUB = SUB + SUBM;
      ZZ = ZZ + 1.;
      goto label2;
    }
  }
  // Expansion coefficients for the first derivative
  double A1 =  1.;
  double A2 =  1./2.;
  double A3 =  1./6.;
  double A4 = -1./30.;
  double A5 =  1./42.;
  double A6 = -1./30.;
  double A7 =  5./66.;
    // Expansion coefficients for the higher derivatives
  if(M>1) {
    for(K2=2;K2<=M;K2++){
      A1 = A1 * (1.*K2-1.);
      A2 = A2 *  1.*K2;
      A3 = A3 * (1.*K2+1.);
      A4 = A4 * (1.*K2+3.);
      A5 = A5 * (1.*K2+5.);
      A6 = A6 * (1.*K2+7.);
      A7 = A7 * (1.*K2+9.);
    }
  }
  dcomplex RZ = 1./ ZZ;
  dcomplex DZ = RZ * RZ;
  // SUB + asympt expansion (Abramowitz, Stengun, 6.4.11)
  return SUB + std::pow(-1.,M+1) * pow(RZ,M) * ( A1 + RZ * (A2 + RZ * (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) );
}

dcomplex psi(dcomplex x) {
  return psi0(x);
}
dcomplex psi(int n, dcomplex x) {
  return psi_n(n,x);
}


////////////////////////////////////////
// Derivatives of Delta(N) = 1/Gamma(N)
//
// GDk means Gamma(N) Delta^{(k)}(N)
void GD(int p, double N, double *lGDk) {
  lGDk[0] = 1;
  for(int k=1; k<=p; k++) {
    double res = 0.;
    for(int i=0; i<k; i++) {
      res -= gsl_sf_choose(k-1,i) * lGDk[i] * gsl_sf_psi_n(k-1-i, N);
    }
    lGDk[k] = res;
  }
}
double GD(int k, double N) {
  double  *lGDk = new double [k+2];
  GD(k+1,N,lGDk);
  double res = lGDk[k];
  delete[] lGDk;
  return res;
}
// complex version
void GD(int p, dcomplex N, dcomplex *lGDk) {
  lGDk[0] = 1;
  for(int k=1; k<=p; k++) {
    dcomplex res = 0.;
    for(int i=0; i<k; i++) {
      res -= gsl_sf_choose(k-1,i) * lGDk[i] * psi_n(k-1-i, N);
    }
    lGDk[k] = res;
  }
}
//////////////////////////////////////////







////////////////////////////////////////
// Derivatives of Gamma(N)
//
// DGk means Delta(N) Gamma^{(k)}(N)
void DG(int p, double N, double *lDGk) {
  lDGk[0] = 1;
  for(int k=1; k<=p; k++) {
    double res = 0.;
    for(int i=0; i<k; i++) {
      res += gsl_sf_choose(k-1,i) * lDGk[i] * gsl_sf_psi_n(k-1-i, N);
    }
    lDGk[k] = res;
  }
}
double DG(int k, double N) {
  double  *lDGk = new double [k+2];
  DG(k+1,N,lDGk);
  double res = lDGk[k];
  delete[] lDGk;
  return res;
}
// complex version
void DG(int p, dcomplex N, dcomplex *lDGk) {
  lDGk[0] = 1;
  for(int k=1; k<=p; k++) {
    dcomplex res = 0.;
    for(int i=0; i<k; i++) {
      res += gsl_sf_choose(k-1,i) * lDGk[i] * psi_n(k-1-i, N);
    }
    lDGk[k] = res;
  }
}
//////////////////////////////////////////






////////////////////////////////////////
// Derivatives of G(N,x) = Gamma(N-x/2) Delta(N+x/2) wrt x
//
// GBk means k-th derivative of G(N,x) in x=0
void GB(int p, double N, double *lGBk) {
  lGBk[0] = 1;
  for(int k=1; k<=p; k++) {
    double res = 0.;
    for(int i=0; i<k; i++) {
      if(i%2 != 0) continue;
      res -= gsl_sf_choose(k-1,i) * lGBk[k-1-i] * gsl_sf_psi_n(i, N) / std::pow(2.,i);
    }
    lGBk[k] = res;
  }
}
double GB(int k, double N) {
  double  *lGBk = new double [k+2];
  GB(k,N,lGBk);
  double res = lGBk[k];
  delete[] lGBk;
  return res;
}
//////////////////////////////////////////
// complex version
void GB(int p, dcomplex N, dcomplex *lGBk) {
  lGBk[0] = 1;
  for(int k=1; k<=p; k++) {
    dcomplex res = 0.;
    for(int i=0; i<k; i++) {
      if(i%2 != 0) continue;
      res -= gsl_sf_choose(k-1,i) * lGBk[k-1-i] * psi_n(i, N) / std::pow(2.,i);
    }
    lGBk[k] = res;
  }
}
//////////////////////////////////////////






////////////////////////////////////////
// Derivatives of CB(x) = Gamma(1+x) Gamma(1-x/2) Delta(1+x/2) wrt x
//
// double CBk[MAX_ORDER+1];
// double CBrec(int k) {
//   if(k == 0) {
//     CBk[0] = 1.;
//     return 1.;
//   }
//   double res = 0.;
//   for(int i=0; i<k; i++) {
//     res += gsl_sf_choose(k-1,i) * CBk[k-1-i] * gsl_sf_psi_n(i, 1.) * (1. - (1.+pow(-1.,i))/pow(2.,i+1) );
//   }
//   CBk[k] = res;
//   return res;
// }
// double CB(int k) {
//   if(k>MAX_ORDER) errmaxorder();
//   for(int i=0; i<=k; i++) {
//     CBrec(i);
//   }
//   return CBk[k];
// }
//////////////////////////////////////////
void CB(int p, double *lCBk) {
  lCBk[0] = 1;
  for(int k=1; k<=p; k++) {
    double res = 0.;
    for(int i=0; i<k; i++) {
      res += gsl_sf_choose(k-1,i) * lCBk[k-1-i] * gsl_sf_psi_n(i, 1.) * (1. - (1.+std::pow(-1.,i))/std::pow(2.,i+1) );
    }
    lCBk[k] = res;
  }
}
double CB(int k) {
  double  *lCBk = new double [k+2];
  CB(k,lCBk);
  double res = lCBk[k];
  delete[] lCBk;
  return res;
}










//////////////////////////////////////////
// Mellin of Dk(z)  = [ log^k (1-z) / (1-z) ]_+
//
double mell_D(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGDk = new double [k+2];
  DG(k+1,1,lDGk);
  GD(k+1,N,lGDk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGDk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}
// derivative wrt N
double mell_D_der(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGDk = new double [k+3];
  DG(k+1,1,lDGk);
  GD(k+2,N,lGDk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * ( lGDk[k+2-i] + gsl_sf_psi(N)*lGDk[k+1-i] );
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}
// complex version
dcomplex mell_D(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGDk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  GD(k+1,N,lGDk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGDk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}
//////////////////////////////////////////


//////////////////////////////////////////
// same as before, but without constant
//
double mell_D_noconst(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGDk = new double [k+2];
  DG(k+1,1,lDGk);
  GD(k+1,N,lGDk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGDk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}
// complex version
dcomplex mell_D_noconst(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGDk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  GD(k+1,N,lGDk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGDk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}
//////////////////////////////////////////






//////////////////////////////////////////
// Mellin of DMk(z) = [ log^k (log(1/z)) / log(1/z) ]_+
//  (large N limit of mell_D, i.e. the Minimal prescription logs)
//
double mell_DM(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  DG(k+1,1,lDGk);
  for(int i=0; i<=k; i++) {  // doesn't reproduce the constant term, as with the MP
    res += gsl_sf_choose(k+1,i) * lDGk[i] * pow(log(1./N),k+1-i);
  }
  delete[] lDGk;
  return res/(k+1.);
}
dcomplex mell_DM(int k, dcomplex N) {
  dcomplex res = 0.;
  double  *lDGk = new double [k+2];
  DG(k+1,1,lDGk);
  for(int i=0; i<=k; i++) {  // doesn't reproduce the constant term, as with the MP
    res += gsl_sf_choose(k+1,i) * lDGk[i] * pow(log(1./N),k+1-i);
  }
  delete[] lDGk;
  return res/(k+1.);
}
//////////////////////////////////////////






//////////////////////////////////////////
// Mellin of DMconstk(z) = [ log^k (log(1/z)) / log(1/z) ]_+ + delta(1-z) Gamma^(k+1) /(k+1)
//  (large N limit of mell_D, i.e. the Minimal prescription logs, including constants)
//
double mell_DMconst(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  DG(k+1,1,lDGk);
  for(int i=0; i<=k+1; i++) {  // reproduces also the constant term, but it is not obtained with the MP
    res += gsl_sf_choose(k+1,i) * lDGk[i] * pow(log(1./N),k+1-i);
  }
  delete[] lDGk;
  return res/(k+1.);
}
dcomplex mell_DMconst(int k, dcomplex N) {
  dcomplex res = 0.;
  double  *lDGk = new double [k+2];
  DG(k+1,1,lDGk);
  for(int i=0; i<=k+1; i++) {  // doesn't reproduce the constant term, as with the MP
    res += gsl_sf_choose(k+1,i) * lDGk[i] * pow(log(1./N),k+1-i);
  }
  delete[] lDGk;
  return res/(k+1.);
}
//////////////////////////////////////////






//////////////////////////////////////////
// Mellin of DBk(z) = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+'
// where prime denotes that the delta term is the same as in D_k(z)
//
double mell_DB(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGBk = new double [k+2];
  DG(k+1,1,lDGk);
  GB(k+1,N,lGBk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////
// complex version
dcomplex mell_DB(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGBk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  GB(k+1,N,lGBk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////


//////////////////////////////////////////
// same as before, but without constant
//
double mell_DB_noconst(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGBk = new double [k+2];
  DG(k+1,1,lDGk);
  GB(k+1,N,lGBk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////
// complex version
dcomplex mell_DB_noconst(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGBk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  GB(k+1,N,lGBk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i];
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////







//////////////////////////////////////////
// Mellin of DBk(z) - DBtildek = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+' - [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+
// i.e. the constant difference between the two possible defs of the distributions
//
double mell_Bdiff(int k) {
  double *lCBk = new double[k+2];
  CB(k+1,lCBk);
  double res = lCBk[k+1]/(k+1.);
  delete[] lCBk;
  return res;
}
//////////////////////////////////////////







//////////////////////////////////////////
// Mellin of Dbark(z) = [ log^k [(1-z)/z] / (1-z) ]_+'
// where prime denotes that the delta term is the same as in D_k(z)
//
double mell_Dbar(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGBk = new double [k+2];
  DG(k+1,1,lDGk);
  DG(k+1,N,lGBk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i] * pow(-1.,k+1-i);
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////
// complex version
dcomplex mell_Dbar(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGBk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  DG(k+1,N,lGBk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i] * pow(-1.,k+1-i);
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////

double mell_Dbar_der(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGDk = new double [k+3];
  DG(k+1,1,lDGk);
  DG(k+2,N,lGDk);
  for(int i=0; i<=k+1; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * pow(-1.,k+1-i) * ( lGDk[k+2-i] - gsl_sf_psi(N)*lGDk[k+1-i] );
  }
  delete[] lDGk;
  delete[] lGDk;
  return res/(k+1.);
}

//////////////////////////////////////////
// same as before, but without constant
//
double mell_Dbar_noconst(int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+2];
  double  *lGBk = new double [k+2];
  DG(k+1,1,lDGk);
  DG(k+1,N,lGBk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i] * pow(-1.,k+1-i);
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////
// complex version
dcomplex mell_Dbar_noconst(int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+2];
  dcomplex *lGBk = new dcomplex[k+2];
  DG(k+1,1,lDGk);
  DG(k+1,N,lGBk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k+1,i) * lDGk[i] * lGBk[k+1-i] * pow(-1.,k+1-i);
  }
  delete[] lDGk;
  delete[] lGBk;
  return res/(k+1.);
}
//////////////////////////////////////////




//////////////////////////////////////////
// Mellin of Lpk(z) = (1-z)^p log^k (1-z)
//
double mell_L(unsigned int p, unsigned int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+1];
  double  *lGDk = new double [k+1];
  DG(k,1+p,lDGk);
  GD(k,N+1+p,lGDk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k,i) * lDGk[i] * lGDk[k-i];
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lGDk;
  return res*gsl_sf_fact(p);
}
// complex version
dcomplex mell_L(unsigned int p, unsigned int k, dcomplex N) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+1];
  dcomplex *lGDk = new dcomplex[k+1];
  DG(k,1+p,lDGk);
  GD(k,N+1+p,lGDk);
  for(int i=0; i<=k; i++) {
    res += gsl_sf_choose(k,i) * lDGk[i] * lGDk[k-i];
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lGDk;
  return res*gsl_sf_fact(p);
}
//////////////////////////////////////////

// hack for soft expansion
dcomplex mell_L_trunc(unsigned int p, unsigned int k, dcomplex N, unsigned int trunc) {
  dcomplex res = 0.;
  double   *lDGk = new double  [k+1];
  dcomplex *lGDk = new dcomplex[k+1];
  DG(k,1+p,lDGk);
  GD(k,N+1+p,lGDk);
  for(int i=0; i<=trunc; i++) {
    res += gsl_sf_choose(k,i) * lDGk[i] * lGDk[k-i];
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lGDk;
  return res*gsl_sf_fact(p);
}
double mell_L_trunc(unsigned int p, unsigned int k, double N, unsigned int trunc) {
  double res = 0.;
  double *lDGk = new double[k+1];
  double *lGDk = new double[k+1];
  DG(k,1+p,lDGk);
  GD(k,N+1+p,lGDk);
  for(int i=0; i<=trunc; i++) {
    res += gsl_sf_choose(k,i) * lDGk[i] * lGDk[k-i];
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lGDk;
  return res*gsl_sf_fact(p);
}



double mell_Lz(unsigned int p, unsigned int k, double N) {
  double res = 0.;
  double  *lDGk = new double [k+1];
  double  *lDGNk = new double [k+1];
  double  *lGDNk = new double [k+1];
  DG(k,1+p,lDGk);
  DG(k,N,lDGNk);
  GD(k,N+1+p,lGDNk);
  for(int i=0; i<=k; i++) {
    double r = 0;
    for(int n=0; n<=k-i; n++) {
      r += gsl_sf_choose(k-i,n) * pow(-1.,n) * lDGNk[n] * lGDNk[k-i-n];
    }
    res += gsl_sf_choose(k,i) * lDGk[i] * r /pow(2.,k-i);
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lDGNk;
  delete[] lGDNk;
  return res*gsl_sf_fact(p);
}
// hack for soft expansion
double mell_Lz_trunc(unsigned int p, unsigned int k, double N, unsigned int trunc) {
  double res = 0.;
  double  *lDGk = new double [k+1];
  double  *lDGNk = new double [k+1];
  double  *lGDNk = new double [k+1];
  DG(k,1+p,lDGk);
  DG(k,N,lDGNk);
  GD(k,N+1+p,lGDNk);
  for(int i=0; i<=trunc; i++) {
    double r = 0;
    for(int n=0; n<=k-i; n++) {
      r += gsl_sf_choose(k-i,n) * pow(-1.,n) * lDGNk[n] * lGDNk[k-i-n];
    }
    res += gsl_sf_choose(k,i) * lDGk[i] * r /pow(2.,k-i);
  }
  for(int i=0; i<=p; i++) {
    res /= N+i;
  }
  delete[] lDGk;
  delete[] lDGNk;
  delete[] lGDNk;
  return res*gsl_sf_fact(p);
}



