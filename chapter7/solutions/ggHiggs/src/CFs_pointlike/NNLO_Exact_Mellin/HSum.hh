#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>
#include <sstream>
#include <string>
#include <functional>
#include "../../math/complex.hh"

using namespace std;

class HSum {	
public:
  HSum(bool verbose=false,bool testinterfun=false, bool testharmsums=false);
  virtual ~HSum();
  dcomplex HS(int i, dcomplex N);
  dcomplex HS(int i, int j, dcomplex N);
  dcomplex HS(int i, int j, int k, dcomplex N);
  dcomplex HS(int i, int j, int k, int m, dcomplex N);
  
  
  
private:
  bool _testinterpolatedfunction;
  bool _testharmonicsums;
  bool _verbose;
  
  void InizializeConst();
  
  //crosscheck function
  
  void TestInterpolatedFunction();
  void TestHarmonicSums(int n);
  dcomplex HS_int(int i, int N);
  dcomplex HS_int(int i, int j, int N);
  dcomplex HS_int(int i, int j, int k, int N);
  dcomplex HS_int(int i, int j, int k, int m, int N);
  
  //Semi-analitics Mellin transform
  dcomplex g1(dcomplex N);
  dcomplex g2(dcomplex N);
  dcomplex g3(dcomplex N);
  dcomplex g4(dcomplex N);
  dcomplex g5(dcomplex N);
  dcomplex g6(dcomplex N);
  dcomplex g7(dcomplex N);
  dcomplex g8(dcomplex N);
  dcomplex g9(dcomplex N);
  dcomplex g10(dcomplex N);
  dcomplex g11(dcomplex N);
  dcomplex g12(dcomplex N);
  dcomplex g13(dcomplex N);
  dcomplex g14(dcomplex N);
  dcomplex g15(dcomplex N);
  dcomplex g16(dcomplex N);
  dcomplex g17(dcomplex N);
  dcomplex g18(dcomplex N);
  dcomplex g19(dcomplex N);
  dcomplex g20(dcomplex N);
  dcomplex g21(dcomplex N);
  dcomplex g22(dcomplex N);
  dcomplex g23(dcomplex N);
  dcomplex g24(dcomplex N);
  dcomplex g25(dcomplex N);
  dcomplex g26(dcomplex N);
  dcomplex g27(dcomplex N);
  dcomplex g28(dcomplex N);
  dcomplex g29(dcomplex N);
  dcomplex g30(dcomplex N);
  dcomplex g31(dcomplex N);
  dcomplex g32(dcomplex N);
  dcomplex g33(dcomplex N);
  dcomplex g34(dcomplex N);
  dcomplex g35(dcomplex N);
  dcomplex g36(dcomplex N);
  dcomplex g37(dcomplex N);
  dcomplex g38(dcomplex N);
  dcomplex g39(dcomplex N);
  
  
  
  
  //Math important Function
  
  
  dcomplex Log(dcomplex N){
  return(log(N));
  }
  
  long double Zeta(int i){
  return gsl_sf_zeta_int(i);
  }

  //Complex Euler Gamma

  dcomplex LogGamma(dcomplex z){
 int g=7;
 double p[9];
 p[0]=0.99999999999980993;
 p[1]=676.5203681218851;
 p[2]=-1259.1392167224028;
 p[3]=771.32342877765313;
 p[4]=-176.61502916214059;
 p[5]=12.507343278686905;
 p[6]=-0.13857109526572012;
 p[7]=9.9843695780195716e-6;
 p[8]=1.5056327351493116e-7;
 dcomplex sum,ris;
 z -=1;
 sum=p[0];
 for(int i=1;i<(g+2);i++)
   sum += p[i]/(z+i);
 ris=0.5*gsl_sf_log(2*M_PI)+(z+0.5)*log(z+g+0.5)-(z+g+0.5)+log(sum);
 return ris;
  }

  dcomplex CGamma(dcomplex z){
  return( exp(LogGamma(z)));
  }

  //Gamma Derivatives of order i complex
  dcomplex PolyGamma(int i, dcomplex z){
  if (i == 0) {
    dcomplex SUB = 0. ;
    dcomplex ZZ = z;
    if(std::abs(imag(ZZ))<10.) { // if too close to the real axis...
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
    return (SUB + log(ZZ) - 0.5 * RZ - DZ/5040. * ( 420.+ DZ * ( - 42. + DZ * (20. - 21. * DZ) ) ));
  }
  else {
     int K1, K2;
     dcomplex SUB = 0. , SUBM;
     dcomplex ZZ = z;
     if(std::abs(imag(ZZ))<10.) { // if too close to the real axis...
      label2:
      SUBM = -1./ZZ;
      for(K1=1; K1<=i; K1++) {
	SUBM = - SUBM * K1 / ZZ;
      }
      if(real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
	SUB = SUB + SUBM;
	ZZ = ZZ + 1.;
	goto label2;
      }
    }
    // Expansion coefficients for the first derivative
    long double A1 =  1.;
    long double A2 =  1./2.;
    long double A3 =  1./6.;
    long double A4 = -1./30.;
    long double A5 =  1./42.;
    long double A6 = -1./30.;
    long double A7 =  5./66.;
    // Expansion coefficients for the higher derivatives
    if(i>1) {
      for(K2=2;K2<=i;K2++){
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
    return (SUB + pow(-1.,i+1) * pow(RZ,i) * ( A1 + RZ * (A2 + RZ * (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) ));
  }
  }
  
  dcomplex B(int i,dcomplex z){
    return(1./pow(2.,(long double)i+1.)*(PolyGamma(i,(1.+z)/2.)-PolyGamma(i,z/2.)));
  }
  
  
  //Important Costants
  long double zeta2;
  long double zeta3;
  long double Li4;
  long double log2;
  long double log2q;
  long double log2c;
  long double zeta2q;
  long double EulerGamma;
  
  
  
  //Vectors of coefficients
  std::vector<long double> a1;
  std::vector<long double> a2;
  std::vector<long double> a3;
  std::vector<long double> b1;
  std::vector<long double> b2;
  std::vector<long double> b3;
  std::vector<long double> c1;
  std::vector<long double> P21;
  std::vector<long double> c2;
  std::vector<long double> P22;
  std::vector<long double> c3;
  std::vector<long double> P23;
  std::vector<long double> P33;
  std::vector<long double> c4;
  std::vector<long double> P24;
  std::vector<long double> P34;
  std::vector<long double> c5;
  std::vector<long double> d5;
  std::vector<long double> q1;
  std::vector<long double> q2;
  std::vector<long double> q3;
  std::vector<long double> q4;
  std::vector<long double> q5;
  std::vector<long double> q6;
  std::vector<long double> q7;
  
  //Table of Harmonic Sums up to weight 4
  //Weight 1
  dcomplex H_1(dcomplex N);
  dcomplex H_m1(dcomplex N);
  //Weight 2
  dcomplex H_2(dcomplex N);
  dcomplex H_m2(dcomplex N);
  dcomplex H_1_1(dcomplex N);
  dcomplex H_1_m1(dcomplex N);
  dcomplex H_m1_1(dcomplex N);
  dcomplex H_m1_m1(dcomplex N);
  //weight 3
  dcomplex H_3(dcomplex N);
  dcomplex H_m3(dcomplex N);
  dcomplex H_m2_m1(dcomplex N);
  dcomplex H_m2_1(dcomplex N);
  dcomplex H_m1_m2(dcomplex N);
  dcomplex H_m1_2(dcomplex N);
  dcomplex H_1_m2(dcomplex N);
  dcomplex H_1_2(dcomplex N);
  dcomplex H_2_1(dcomplex N);
  dcomplex H_2_m1(dcomplex N);
  dcomplex H_1_1_1(dcomplex N);
  dcomplex H_1_1_m1(dcomplex N);
  dcomplex H_1_m1_1(dcomplex N);
  dcomplex H_1_m1_m1(dcomplex N);
  dcomplex H_m1_1_1(dcomplex N);
  dcomplex H_m1_1_m1(dcomplex N);
  dcomplex H_m1_m1_1(dcomplex N);
  dcomplex H_m1_m1_m1(dcomplex N);
  //Weight 4
  dcomplex H_4(dcomplex N);
  dcomplex H_m4(dcomplex N);
  dcomplex H_m3_m1(dcomplex N);
  dcomplex H_m3_1(dcomplex N);
  dcomplex H_m1_m3(dcomplex N);
  dcomplex H_m1_3(dcomplex N);
  dcomplex H_1_m3(dcomplex N);
  dcomplex H_1_3(dcomplex N);
  dcomplex H_3_1(dcomplex N);
  dcomplex H_3_m1(dcomplex N);
  dcomplex H_m2_m2(dcomplex N);
  dcomplex H_m2_2(dcomplex N);
  dcomplex H_2_m2(dcomplex N);
  dcomplex H_2_2(dcomplex N);
  dcomplex H_m2_1_1(dcomplex N);
  dcomplex H_m2_1_m1(dcomplex N);
  dcomplex H_m2_m1_1(dcomplex N);
  dcomplex H_m2_m1_m1(dcomplex N);
  dcomplex H_m1_1_2(dcomplex N);
  dcomplex H_m1_1_m2(dcomplex N);
  dcomplex H_m1_m1_2(dcomplex N);
  dcomplex H_m1_m1_m2(dcomplex N);
  dcomplex H_m1_2_1(dcomplex N);
  dcomplex H_m1_2_m1(dcomplex N);
  dcomplex H_m1_m2_1(dcomplex N);
  dcomplex H_m1_m2_m1(dcomplex N);
  dcomplex H_1_1_2(dcomplex N);
  dcomplex H_1_1_m2(dcomplex N);
  dcomplex H_1_m1_2(dcomplex N);
  dcomplex H_1_m1_m2(dcomplex N);
  dcomplex H_1_2_1(dcomplex N);
  dcomplex H_1_2_m1(dcomplex N);
  dcomplex H_1_m2_1(dcomplex N);
  dcomplex H_1_m2_m1(dcomplex N);
  dcomplex H_2_1_1(dcomplex N);
  dcomplex H_2_1_m1(dcomplex N);
  dcomplex H_2_m1_1(dcomplex N);
  dcomplex H_2_m1_m1(dcomplex N);
  dcomplex H_1_1_1_1(dcomplex N);
  dcomplex H_1_1_1_m1(dcomplex N);
  dcomplex H_1_1_m1_1(dcomplex N);
  dcomplex H_1_1_m1_m1(dcomplex N);
  dcomplex H_1_m1_1_1(dcomplex N);
  dcomplex H_1_m1_1_m1(dcomplex N);
  dcomplex H_1_m1_m1_1(dcomplex N);
  dcomplex H_1_m1_m1_m1(dcomplex N);
  dcomplex H_m1_1_1_1(dcomplex N);
  dcomplex H_m1_1_1_m1(dcomplex N);
  dcomplex H_m1_1_m1_1(dcomplex N);
  dcomplex H_m1_1_m1_m1(dcomplex N);
  dcomplex H_m1_m1_1_1(dcomplex N);
  dcomplex H_m1_m1_1_m1(dcomplex N);
  dcomplex H_m1_m1_m1_1(dcomplex N);
  dcomplex H_m1_m1_m1_m1(dcomplex N);
  
 
  
  
};

