#pragma once
#include<string>

namespace ggHiggs {

  extern int pdf_init(std::string SET_NAME, int SUBSET);
  extern double alphasPDF(double mu);
  extern double getQuarkMass(int flavour);

  extern double lumi_gg (double x1, double x2, double mu);
  extern double lumi_qg (double x1, double x2, double mu);
  extern double lumi_qqb(double x1, double x2, double mu);
  extern double lumi_qq (double x1, double x2, double mu);
  extern double lumi_qqp(double x1, double x2, double mu);
  extern double lumi_SS (double x1, double x2, double mu);
  extern double lumi_VV (double x1, double x2, double mu);
  extern double lumi_TT (double x1, double x2, double mu);

  inline double lumi(double lum(double,double,double), double t, double x, double mu) {
    double x1 = 1.-(1.-x)*t;
    double x2 = x/x1;
    return lum(x1, x2, mu) /x1 * (1-x);
  }

  extern double Lumij(double tau, double mu, int LumChannel);
  extern double Lum(double tau, double mu);
  extern double ChebRecLum(double *c, int N, double z, double beta, int gamma);

};
