#pragma once

#include "HSum.hh"
#include "../../math/complex.hh"
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

class MellinFunc{
public:
  MellinFunc();
  dcomplex Li3zplus(dcomplex N); // Li3(z)/(1+z)
  dcomplex S12z2plus(dcomplex N); // S12(z^2)/(1+z)
  dcomplex S12mzplus(dcomplex N); // S12(-z)/(1+z)
  dcomplex S12zplus(dcomplex N); // S12(z)/(1+z)
  dcomplex Li3mzplus(dcomplex N); //Li3(-z)/(1+z)
  dcomplex Li2zLogminus(dcomplex N); //Li2(z)Log(1-z)
  dcomplex Li2zLogplus(dcomplex N); //Li2(z)Log(1+z)
  dcomplex Li2mzLogplus(dcomplex N); //Li2(-z)Log(1+z)
  dcomplex Li2mzLogz(dcomplex N); //Li2(-z)Log(z)
  dcomplex Li2zLogzplusminus(dcomplex N); //Li2(z)Log(z)/((1+z)(1-z))
  dcomplex Li2mzLogzplus(dcomplex N); //Li2(-z)Log(z)/(1+z)
  dcomplex Li2mzLogminus2plus(dcomplex N); //Li2(-z)Log[1-z]^2/(1+z)
  dcomplex Li2mzLogplusplus(dcomplex N); //Li2(-z)Log[1+z]/(1+z)
  dcomplex Logz2Logminusminus(dcomplex N); //Log(z)^2Log(1-z)/(1-z)
  dcomplex LogzLogminus2minus(dcomplex N); //Log(z)Log(1-z)^2/(1-z)
  dcomplex Liz2Logminus(dcomplex N); //Li2(z)Log(1-z)
  dcomplex Logz3minusplus(dcomplex N); // Log(z)^3 /((1-z)(1+z))
  dcomplex Logplusplus(dcomplex N); //Log(1+z)/(1+z)
  dcomplex Li2zLogplusplus(dcomplex N); //Li2(z)Log(1+z)/(1+z)
  dcomplex Logz2Logplusplus(dcomplex N); //Log(1+z)Log(z)^2/(1+z)
  dcomplex LogzLogplus(dcomplex N); //Log(z)Log(1+z)
  dcomplex Logz2Logplus(dcomplex N); //Log(z)^2Log(1+z)
  dcomplex LogzLogplus2(dcomplex N); //Log(z)Log(1+z)^2
  dcomplex LogzLogplus2plus(dcomplex N); //Log(z)Log(1+z)^2/(1+z)
  dcomplex Li2z(dcomplex N); //Li2(z)
  dcomplex Li2mz(dcomplex N); //Li2(-z)
  dcomplex S12z(dcomplex N); //S12(z)
  dcomplex S12mz(dcomplex N); //S12(mz)
  dcomplex S12z2(dcomplex N); //S12(z^2)
  dcomplex Li3z(dcomplex N); //Li3(z)
  dcomplex Li3mz(dcomplex N); //Li3(mz)
  dcomplex Li2zLogz(dcomplex N); //Li2(z)Log(z)
  dcomplex Logminus(dcomplex N); //Log(1-z)
  dcomplex Logplus(dcomplex N); //Log(1+z)
  dcomplex Logz(dcomplex N); //Log(z)
  dcomplex Logz2(dcomplex N); //Log(z)^2
  dcomplex Logz3(dcomplex N); //Log(z)^3
  dcomplex LogzLogminus(dcomplex N); //Log(z)Log(1-z)
  dcomplex LogzLogminus2(dcomplex N); //Log(z)Log(1-z)^2
  dcomplex Logz2Logminus(dcomplex N); //Log(z)^2Log(1-z)
  dcomplex Logz2minus(dcomplex N); //log(z)^2/(1-z)
  dcomplex Logminus2(dcomplex N); //Log(1-z)^2
  dcomplex Logminus3(dcomplex N); //Log(1-z)^3
  dcomplex LogzLogminusminus(dcomplex N); //Log(z)Log(1-z)/(1-z)
  dcomplex Logzminus(dcomplex N); //log(z)/(1-z)
  dcomplex Logzminusplus(dcomplex N); //log(z)/((1-z)(1+z))
  dcomplex plus(dcomplex N); //1/(1+z)
  dcomplex Li2minusminus(dcomplex N); //Li2(1-z)/(1-z)
  dcomplex S12zregminus(dcomplex N); //(S12(z)-Zeta(3))/(1-z)
  dcomplex Li2zregminus(dcomplex N); //(Li2(z)-Zeta(2))/(1-z)
  dcomplex Li3zregminus(dcomplex N); //(Li3(z)-Zeta(3))/(1-z)
  dcomplex Li2zregLogminusminus(dcomplex N); //(Li2(z)-Zeta(2))Log(1-z)/(1-z)
  dcomplex D0(dcomplex N); //[1/(1-z)]_+
  dcomplex D1(dcomplex N); //[Log(1-z)/(1-z)]_+
  dcomplex D2(dcomplex N); //[Log(1-z)^2/(1-z)]_+
  dcomplex D3(dcomplex N); //[Log(1-z)^3/(1-z)]_+
  
  
  
  
private:
  HSum H;
  long double zeta2;
  long double zeta3;
  long double Li4;
  long double log2;
  long double log2q;
  long double log2c;
  long double zeta2q;
  long double EulerGamma;
  
};

dcomplex sigma_NNLO_gg_N(dcomplex );
dcomplex sigma_NNLO_qg_N(dcomplex );
dcomplex sigma_NNLO_qqbar_N(dcomplex );
dcomplex sigma_NNLO_qq_N(dcomplex );
dcomplex sigma_NNLO_qqprime_N(dcomplex );


