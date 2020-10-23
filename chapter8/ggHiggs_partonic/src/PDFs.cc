#include "LHAPDF/LHAPDF.h"  // interface with LHAPDF
#include "PDFs.hh"
#include "parameters.hh"
#include "math/integration.hh"
#include <gsl/gsl_sf.h> 
#include <gsl/gsl_integration.h>

using namespace std;

namespace ggHiggs {

  const int _nf = 5;  // this is the max value of nf

#if LHAPDF_MAJOR_VERSION == 6
  LHAPDF::PDF* lhpdf;
#endif

  double xfx(double x, double Q, int parton) {
    //if(excludePDF[6+parton]) return 0;
#if LHAPDF_MAJOR_VERSION == 6
    //if(parton == 0) parton = 21;
    return lhpdf->xfxQ(parton, x, Q);
#else
    return LHAPDF::xfx(x,Q,parton);
#endif
  }



  void xfx(double x, double Q, double *pdfs) {
#if LHAPDF_MAJOR_VERSION == 6
    vector<double> rtn;
    lhpdf->xfxQ(x, Q, rtn);
    for(int i=0; i<=12; i++) {
      //pdfs[i] = ( excludePDF[i] ? 0 : rtn[i] );
      pdfs[i] = rtn[i];
    }
#else
    LHAPDF::xfx(x,Q,pdfs);
    //for(int i=0; i<=12; i++) if(excludePDF[i]) pdfs[i] = 0;
#endif
    return;
  }


/////////////////////////////
//  LUMINOSITY
double lumi_gg(double x1, double x2, double mu) {
  //if(x1<1.e-7) x1=1.e-7;
  //if(x2<1.e-7) x2=1.e-7;
  return xfx(x1,mu,0) * xfx(x2,mu,0) /x1/x2;
}
double lumi_qg(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double res = 0;
  for(int i=1; i<=_nf; i++) {
    res += (xfx1[6+i] + xfx1[6-i]) * xfx2[6+0];
    res += (xfx2[6+i] + xfx2[6-i]) * xfx1[6+0];
  }
  return res/x1/x2;
}
double lumi_qqb(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double res = 0;
  for(int i=-_nf; i<=_nf; i++) {
    if(i==0) continue;
    res += xfx1[6+i] * xfx2[6-i];
  }
  return res/x1/x2;
}
double lumi_qq(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double res = 0;
  for(int i=-_nf; i<=_nf; i++) {
    if(i==0) continue;
    res += xfx1[6+i] * xfx2[6+i];
  }
  return res/x1/x2;
}
double lumi_qqp(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double res = 0;
  for(int i=-_nf; i<=_nf; i++) {
    for(int j=-_nf; j<=_nf; j++) {
      if(i!=j && i!=-j && i!=0 && j!=0) res += xfx1[6+i] * xfx2[6+j];
    }
  }
  return res/x1/x2;
}
double lumi_SS(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double S1=0, S2=0;
  for(int i=1; i<=_nf; i++) {
    S1 += xfx1[6+i]+xfx1[6-i];
    S2 += xfx2[6+i]+xfx2[6-i];
  }
  return S1*S2/x1/x2;
}
// this is = Vu^2 + Vd^2 + Vs^2 + ...
double lumi_VV(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double VV=0, V1=0, V2=0;
  for(int i=1; i<=_nf; i++) {
    V1 = xfx1[6+i]-xfx1[6-i];
    V2 = xfx2[6+i]-xfx2[6-i];
    VV += V1*V2;
  }
  return VV/x1/x2;
}
/*
// this is = V^2/nf + 1/2 [ V3^2 + V8^2/3 + V15^2/6 + V24^2/10 (+ V35^2/15) ]
double lumi_VV(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double Vi1=0, Vi2=0, VV=0, V1=0, V2=0;
  for(int i=1; i<=_nf; i++) {
    V1 += xfx1[6+i]-xfx1[6-i];
    V2 += xfx2[6+i]-xfx2[6-i];
  }
  for(int i=2; i<=_nf; i++) {
    Vi1=0; Vi2=0;
    for(int j=1; j<=i-1; j++) {
      Vi1 += xfx1[6+j]-xfx1[6-j];
      Vi2 += xfx2[6+j]-xfx2[6-j];
    }
    Vi1 -= (i-1)*(xfx1[6+i]-xfx1[6-i]);
    Vi2 -= (i-1)*(xfx2[6+i]-xfx2[6-i]);
    VV += Vi1*Vi2/(i-1.)/double(i);
  }
  return (VV + V1*V2/double(_nf))/x1/x2;
}
*/
// this is = 1/2 [ T3^2 + T8^2/3 + T15^2/6 + T24^2/10 (+ T35^2/15) ]
double lumi_TT(double x1, double x2, double mu) {
  double xfx1[13], xfx2[13];
  xfx(x1, mu, xfx1);
  xfx(x2, mu, xfx2);
  double T1=0, T2=0, TT=0;
  for(int i=2; i<=_nf; i++) {
    T1=0; T2=0;
    for(int j=1; j<=i-1; j++) {
      T1 += xfx1[6+j]+xfx1[6-j];
      T2 += xfx2[6+j]+xfx2[6-j];
    }
    T1 -= (i-1)*(xfx1[6+i]+xfx1[6-i]);
    T2 -= (i-1)*(xfx2[6+i]+xfx2[6-i]);
    TT += T1*T2/(i-1.)/double(i);
  }
  return TT/x1/x2;
}
/////////////////







////////////////////////////////////////
// integrating delta term
  struct int1dp { double tau, mu; int channel; };
  double Lumij_integrand(double x, void *par) {
    int1dp p = *(int1dp *) par;
    double tau = p.tau;
    double mu  = p.mu;
    int LumChannel = p.channel;
    double res = 0;
    if(LumChannel==0) res = lumi_gg (tau/x,x,mu);
    if(LumChannel==1) res = lumi_qg (tau/x,x,mu);
    if(LumChannel==2) res = lumi_SS (tau/x,x,mu); // Singlet-Singlet
    if(LumChannel==3) res = lumi_TT (tau/x,x,mu);
    if(LumChannel==4) res = lumi_VV (tau/x,x,mu);
    return res/x;
  }
  double Lumij(double tau, double mu, int LumChannel) {
    int1dp p = { tau, mu, LumChannel };
    double prec = 5.e-5;
    gsl_error_handler_t *old_handler = gsl_set_error_handler (NULL);
    gsl_set_error_handler_off();
    double res = integration::gauss(Lumij_integrand, tau, 1., prec, NULL, &p);
    gsl_set_error_handler (old_handler);
    return res;
  }
  double Lum(double tau, double mu) {
    return Lumij(tau,mu,0);
  }
  double ChebRecLum(double *c, int N, double z, double beta, int gamma) {
    double value = 0;
    for(int i=0; i<=N; i++) {
      for(int j=0; j<=gamma; j++) {
	value += c[i] * pow(z,j-beta)* pow(-1.,j) * gsl_sf_choose(gamma,j) *pow(-log(z),i)/gsl_sf_fact(i);
      }
    }
    return value;
  }


  // interface to LHAPDF
  int pdf_init(string PDFset, int member) {
    double asmz;
#if LHAPDF_MAJOR_VERSION == 6
    lhpdf = LHAPDF::mkPDF(PDFset, member);
    LHAPDF::setVerbosity(0);
    asmz = lhpdf->alphasQ(mZ);
#else
    PDFset = PDFset + ".LHgrid";
    LHAPDF::setVerbosity(LHAPDF::SILENT);
    LHAPDF::initPDFSet(PDFset, member);
    asmz = LHAPDF::alphasPDF(mZ);
#endif
    cout << ">>> PDF set:  " << PDFset << " ,  member: " << member
	 << " ,  with alpha_s(mZ) = " << asmz << endl;
    return 0;
  }

  // alphas from LHAPDF
  double alphasPDF(double mu) {
#if LHAPDF_MAJOR_VERSION == 6
    return lhpdf->alphasQ(mu);
#else
    return LHAPDF::alphasPDF(mu);
#endif
  }

  // get Quark masses from LHAPDF
  double getQuarkMass(int flavour) {
    flavour = fabs(flavour);
#if LHAPDF_MAJOR_VERSION == 6
#if LHAPDF_VERSION_CODE >= 60106
    return lhpdf->quarkMass(flavour);
#else
    string flav[] = {"","MDown","MUp","MStrange","MCharm","MBottom","MTop"};
    return lhpdf->info().get_entry_as<double>(flav[flavour]);
#endif
#else
    return LHAPDF::getQMass(flavour);
#endif
  }

};
