#ifndef __ggHiggs__
#define __ggHiggs__

#include "small-x/finite-mt_small-x_coeffs.hh"
#include "PDFs.hh"

using namespace std;


namespace ggHiggs {

  double GetZmass();
  extern double mtop, mtmt;
  extern double mbot; // negative mass means no bottom running in the loop
  extern double mcha; // negative mass means no charm running in the loop
  extern int    order;
  //extern int PDForder;
  extern bool finite_mt;
  extern bool __all_channels__;
  extern SmallxCoeffs *sxc;
  extern double _muF_mu_ratio;
  extern double _muR_muF_ratio;
  extern bool MSbarMass;
  extern bool MSbar;
  extern bool quiet;
  extern double INTEGRAL_PRECISION;
  extern bool isPseudoscalar;


  // actual alphas routine will be provided by the user
  extern double alphas(double mu);

  double alphas_ihixs(double mu, int porder, double as0, double mu0); // alpha_s running from ihixs
  void init_consts(double mH);
  //void SetScales(double mH) {};  // fake - kept it for compatibility issues
  void SetScaleDependence(bool exact_muR=true, bool exact_muF=true);
  void SetSmallxMode(int mode=0, int logstozero=0); // at fixed order
  void SetSmallxResummation(bool active=true, string path="");
  void SetSmallxResummationMode(bool use_LLp=false, int variation=0); // default false
  double Lum(double tau, double muF);
  double xs_LO_prefactor(double mH);
  double xs_LO_prefactor_EFT();
  //double EFTrescaling_deriv(double mH, double mt, double mb, double mc, int id, double &second_deriv, double &third_deriv);
  double EFT_rescaling_MSbar_NLO(double mH);

  void compute(double mH, double sqrts, double LO, double NLO, double NNLO, double *rexact);
  //
  double computeSxRes(double mH, double sqrts, double LO, int ord, double inputalphas=-1);
  //
  // interface to LHAPDF
  int pdf_init(string SET_NAME, int SUBSET);

  double deltaEW(double mH, double mt);

  double FO_3_smallx(double z, double mu, int truncation=37);
  void SetN3LOmode(int mode=2, int truncationOrder=37);

  enum CPcharge { CPeven=0, CPodd=1};
  void SetCPcharge(CPcharge cp);

};


#endif
