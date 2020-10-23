#include <iomanip>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
//#include <string>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sys/stat.h> 
#include <sys/time.h> 
#include <gsl/gsl_sf.h> 
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <getopt.h> 

#include <cuba.h>

#include "math/complex.hh"
#include "math/mellin_Dk.hh"
#include "math/integration.hh"
#include "math/special_functions.hh"
#include "math/chebyshev.hh"
#include "CFs_pointlike/CFs_pointlike.hh"
#include "NLO_finite-mt/bdv.hh"
#include "NNLO_finite-mt/nnlo_finite-mt.hh"
#include "NNLO_finite-mt/nlo_finite-mt.hh"
#include "ggHiggs.hh"
#include "msbar_mass.hh"
#include "parameters.hh"

#ifdef withHELL
#include <hell-x.hh>
#endif

using namespace std;




namespace ggHiggs {

#include "deltaEW.cc"

  const int _nf = 5;
  double _alphas;

  // PRIVATE flags and variables ( don't change here !!! )
  double mtop = 172.5;  // Pole mass
  double mtmt = 162.7;  // MSbar mass
  double mbot = -1; // negative mass means no bottom running in the loop
  double mcha = -1; // negative mass means no charm running in the loop
  vector<double> qmasses;
  dcomplex _mH2_;
  int order = 1;
  //int PDForder = 1;
  bool finite_mt = false;
  SmallxCoeffs *sxc;
  double _muF_mu_ratio  = 1.;
  double _muR_muF_ratio = 1.;
  double muFpdf, mut;
  bool MSbarMass = false;
  bool MSbar = false;
  bool quiet = true;
  double INTEGRAL_PRECISION = 15.e-3;
  double _integration_v_;

  bool _gg_channel_ = true;
  bool _qg_channel_ = true;
  bool _qq_channel_ = true;

  bool newevolmuFmuR = false;

  bool isPseudoscalar = false;

  void SetCPcharge(CPcharge cp) { isPseudoscalar = (cp==CPodd); }




  // Small-x resummation
  bool smallxRes = false;
  bool computeSmallxRes = false;
  bool useLLp = false;
  int  smallxRCvariation = 0;
  void SetSmallxResummationMode(bool use_LLp, int variation) {
    useLLp = use_LLp;
    smallxRCvariation = variation;
  }
  //
#ifdef withHELL
  HELLx::HELLxnf *sxD, *sxDv;
  HELLx::Order SXmatch[4] = { HELLx::LO, HELLx::NLO, HELLx::NNLO, HELLx::N3LO };
  double deltaCggH(double as, double z, double muFrat, int match) {
    sxDv->SetRCvar(smallxRCvariation);
    sxD ->SetRCvar(smallxRCvariation);
    double res = 0;
    if(useLLp) res = sxDv->deltaCggH(as, z, muFrat, SXmatch[match]);
    else       res = sxD ->deltaCggH(as, z, muFrat, SXmatch[match]);
    return res;
  }
  double deltaCggHaux(double as, double z, double muFrat, int match) {
    sxDv->SetRCvar(smallxRCvariation);
    sxD ->SetRCvar(smallxRCvariation);
    double res = 0;
    if(useLLp) res = sxDv->deltaCggHaux(as, z, muFrat, SXmatch[match]);
    else       res = sxD ->deltaCggHaux(as, z, muFrat, SXmatch[match]);
    return res;
  }
  void SetSmallxResummation(bool active, string path) {
    smallxRes = active;
    if(smallxRes) {
      if(sxD ) delete sxD;
      if(sxDv) delete sxDv;
      if(path=="") path = HELLdata;  // HELLdata is passed as a gcc macro
      sxD  = new HELLx::HELLxnf(_nf, HELLx::NLL, path);
      sxDv = new HELLx::HELLxnf(_nf, HELLx::NLL, path);
      sxD ->SetLLpMode(false);
      sxDv->SetLLpMode(true);
      //sxD->SetRCvar();
      sxD ->SetQuietMode();
      sxDv->SetQuietMode();
    }
  }
#else
  double deltaCggH(double as, double z, double muFrat, int match) { return 0; }
  double deltaCggHaux(double as, double z, double muFrat, int match) { return 0; }
  void SetSmallxResummation(bool active, string path) {
    cout << "\033[0;31m  ERROR: ggHiggs was compiled without HELL-x. Small-x resummation not available.\033[0m  " << endl;
  }
#endif
  //
  double SxRES_gg(double as, double z, double muFrat=1, int match=1) {
    double res = deltaCggH(as, z, muFrat, match);
    return res;
  }
  double SxRES_qg(double as, double z, double muFrat=1, int match=1) {
    double res = CF/CA * ( deltaCggH(as, z, muFrat, match) + deltaCggHaux(as, z, muFrat, match) );
    return res;
  }
  double SxRES_qq(double as, double z, double muFrat=1, int match=1) {
    double res = pow(CF/CA,2) * ( deltaCggH(as, z, muFrat, match) + 2*deltaCggHaux(as, z, muFrat, match) );
    return res;
  }



  // small-x at fixed order
  int SXmode = 0;  /// 0=damping; 1=subtraction
  int SXlogs = 0;  /// 0=include all powers of logs at small x; 1=exclude log^0; 2=eclude also log^1 (not a good idea, since at NNLO that's the LL)
  void SetSmallxMode(int mode, int logstozero) {
    if(mode > 1 || mode < -1) {
      cout << "Error: wrong small-x mode" << endl;
      exit(97);
    }
    SXmode = mode;
    SXlogs = logstozero;
  }


  // expansion of log^2(1/z)/z in powers of (z-1) from 0 to 37
  double log2z_softexp[38] = {0, 0, 1., -2., 2.916666667, -3.75, 4.511111111, -5.211111111, 5.859325397, -6.463293651, 7.029087302, -7.561626984, 8.064939875, -8.542356902, 8.996661725, -9.430203368, 9.844981992, -10.24271481, 10.62488732, -10.99279343, 11.34756740, -11.69020927, 12.02160551, -12.34254580, 12.65373676, -12.95581341, 13.24934866, -13.53486123, 13.81282242, -14.08366181, 14.34777206, -14.60551316, 14.85721599, -15.10318539, 15.34370294, -15.57902922, 15.80940597, -16.03505782 };
  double log1z_softexp[38] = {0, -1., 1.5, -1.833333333, 2.083333333, -2.283333333, 2.45, -2.592857143, 2.717857143, -2.828968254, 2.928968254, -3.019877345, 3.103210678, -3.180133755, 3.251562327, -3.318228993, 3.380728993, -3.439552523, 3.495108078, -3.547739657, 3.597739657, -3.645358705, 3.690813250, -3.734291511, 3.775958178, -3.815958178, 3.854419716, -3.891456753, 3.927171039, -3.961653798, 3.994987131, -4.027245195, 4.058495195, -4.088798226, 4.118209990, -4.146781419, 4.174559197, -4.201586224};
  //
  unsigned long int factorials[] = {1,1,2,6,24,120,720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000 };
  unsigned long int factorial(unsigned int k) {
    if(k<16) return factorials[k];
    else return k*factorial(k-1);
    //double res = 1;
    //for(int j=1; j<=k; j++) res *= j;
    //return res;
  }


  // alpha_s running from ihixs
  double alphas_ihixs(double mu, int porder, double as0, double mu0=mZ) {
    // Nf = 5 **************** <---------
    double bb0 = 1.916666667;
    double bb1 = 1.260869565;
    double bb2 = 1.474788647;
    double bb3 = 9.835916431;
    //
    double as;
    as0 /= M_PI;
    //double as0 = asZ/M_PI;
    //double as0=asZ/3.1415926;
    //double Mz=91.187;
    //
    int N = 100000;
    //double step = 2.*log(mu/mZ)/N, incr;
    //double step = 2.*log(mu/Mz)/N, incr;
    double step = 2.*log(mu/mu0)/N, incr;
    for (int i=0; i<N; i++) {
      if (porder >= 3)
	incr = step*bb0*as0*as0*(1+bb1*as0 + bb2*as0*as0 + bb3*as0*as0*as0 );
      else if (porder == 2)
	incr = step*bb0*as0*as0*(1+bb1*as0 + bb2*as0*as0 );
      else if (porder == 1)
	incr = step*bb0*as0*as0*(1+bb1*as0);
      else if (porder == 0)
	incr = step*bb0*as0*as0;
      as = as0 - incr;
      as0 = as;
    }
    return as*M_PI;
  }







  double b0 = beta0(_nf);
  double b1 = beta1(_nf);
  double b2 = beta2(_nf);


  double e01=CA, e00, e11, e22;
  void initGammaCoeffs(bool gamma_NLL, int nf) {
    if(gamma_NLL) {
      e00 = (-11.*CA+2.*nf*(2.*CF/CA-1.))/12.;
      e11 = nf*(26.*CF-23.*CA)/36.;
      e22 = -CA*( (18.*ZETA2-71.)*(CA-2*CF)*nf + CA*CA*(54*ZETA3+99.*ZETA2-395) ) /108.;
    }
    else e00=e11=e22=0;
  }



double Dk(int k, double z) {
  return pow(log(1-z),k)/(1.-z);
}







// perturbative expansions of constants and logs
double _cd1, _cd2, _cd3;
void init_delta_const(int nf) {
  double lf = -2.*log(_muF_mu_ratio);
  _cd1 = 2.*CA * ZETA2;
  _cd2 = ( 837./16. + 67./2.*ZETA2 - 9./20.*ZETA2*ZETA2 - 165./4.*ZETA3
	   + nf  * ( 5./6.*ZETA3 - 5./3.*ZETA2 - 247./36. )
	   -18.*ZETA2*lf*lf + (-27./2. + 171./2.*ZETA3 + 11./6.*nf - (33.-2.*nf)/2.*ZETA2)*lf );
  _cd3 = 1740.02 - 131.928*_nf + 1.75719*_nf*_nf // new coefficient
    - ((1201./576. +5./54.*6*ZETA2 -5./18.*ZETA3)*nf*nf
       + (-3./8.*6*ZETA2 +395./6.*ZETA3 -29807./576. - 5./96.*pow(6*ZETA2,2.))*nf
       + 18217./64. +1089./2.*3*ZETA3*ZETA2 -5049./2.*ZETA5 -11.*(33-2*nf)/12. +6.*(2857./128.-5033./1152.*nf+325./3456.*nf*nf)
       - 46.*2*ZETA2 +5.*6*ZETA2*(51./8.-19./24.*nf) -9453./8.*ZETA3 +55./64.*pow(6*ZETA2,2.)) * lf
    + ((-2./9. + ZETA2/6.)*nf*nf + (81./4.*ZETA3 + 17./3. -2./3.*(51./8.-19./24.*nf) +9./2.*ZETA2)*nf
       - 2673./8.*ZETA3 - 27./10.*pow(6*ZETA2,2.) - 33. - 415./16.*6*ZETA2 + 11.*(51./8.-19./24.*nf)) * lf*lf
    + (72.*ZETA3 +3./2.*(33-2*nf)*ZETA2) * lf*lf*lf;
}
double _W1, _W2, _W3;
double _W1t, _W2t, _W3t;
void init_Wilson_pseudoscalar(double mH) {
  double Ltop = 2*log(_muF_mu_ratio*mH/mtop);
  double CJ2 = 0.;
  double d1c = CA*(2.+2.*ZETA2);
  double d2c = CA*CA*(494./3+1112./9*ZETA2-4./5.*pow(ZETA2,2)-220./3*ZETA3)
    +   CA*_nf*(-82./3-80./9*ZETA2-8./3*ZETA3)
    +   CF*_nf*(-160./3+12*Ltop+16.*ZETA3);
  d2c /= 16.;
  double d3c = _nf*CJ2*(-4.)
    +CF*_nf*_nf*(1498./9-40./9*ZETA2-32./45*pow(ZETA2,2)-224./3*ZETA3)
    +CF*CF*_nf*(457./3+208.*ZETA3-320.*ZETA5)
    +CA*CA*_nf*(-113366./81-10888./81*ZETA2+21032./135*pow(ZETA2,2)+8840./27*ZETA3-2000./3*ZETA2*ZETA3+6952./9*ZETA5)
    +CA*CA*CA*(114568./27+137756./81*ZETA2-61892./135*pow(ZETA2,2)-64096./105*pow(ZETA2,3)-3932*ZETA3+7832./3*ZETA2*ZETA3+13216./3*pow(ZETA3,2)-30316./9*ZETA5)
    +CA*_nf*_nf*(6914./81-1696./81*ZETA2-608./45*pow(ZETA2,2)+688./27*ZETA3)
    +CA*CF*_nf*(-1797.-4160./9*ZETA2+176./45*pow(ZETA2,2)+1856./3*ZETA3+192.*ZETA2*ZETA3+160.*ZETA5+96.*Ltop*(1+ZETA2));
  d3c /= 64.;
  _W1 = d1c - _cd1;
  _W2 = d2c - _cd2 - _cd1*_W1;
  _W3 = d3c - _cd3 - _cd2*_W1 - _cd1*_W2;
}
void init_Wilson(double mH) {
  if(isPseudoscalar) init_Wilson_pseudoscalar(mH);
  else {
    double lHt = 2.*log(_muF_mu_ratio*mH/mtop);
    _W1 = 11./2.;
    _W2 = 121./16. + 2777./144. + 19./8.*lHt + (-67./48. + 2./3.*lHt)*_nf;
    // this was MSbar
    //_W3 = 40.6488 + 40.1458*lHt + 10.941*lHt*lHt; // write full expression .... done.
    //_W3 = (-1792967./20736. +  897943./4608.*ZETA3 + 5347./288.*lHt + 209./32.*lHt*lHt
    //     + _nf * (493./10368. + 209./54.*lHt + 23./16.*lHt*lHt - 110779./6912.*ZETA3)
    //	 + _nf*_nf * (-6865./15552. + 77./864.*lHt - lHt*lHt/9. ) );
    // this is OS
    _W3 = (-1661639./20736. +  897943./4608.*ZETA3 + 6715./288.*lHt + 209./32.*lHt*lHt
	   + _nf * (18925./10368. + 281./54.*lHt + 23./16.*lHt*lHt - 110779./6912.*ZETA3)
	   + _nf*_nf * (-6865./15552. + 77./864.*lHt - lHt*lHt/9. ) );
    if(MSbarMass) _W3 += (19./8. + 2./3.*_nf)*(-8./3.-4*log(mut/mtmt));  // conversion to MSbar mass
    //cout << _W1/pow(M_PI,1) << endl;
    //cout << _W2/pow(M_PI,2) << endl;
    //cout << _W3/pow(M_PI,3) << endl;
  }
}
void init_Wilson_mt(double mH) {
  _W1t = G2l(mH,mtop,mbot,mcha);
  _W2t = gg_nnlo_mt_delta(mH, mtop, _muF_mu_ratio) - _cd2 - _cd1*_W1t;
  _W3t = _W3; // pointlike value
  if(mH > 2*mtop) _W2t = _W2; // pointlike value
}











void init_consts(double mH) {
  if(MSbarMass) mtop = runmass(mtmt, alphas(mtmt)/M_PI, alphas(mut)/M_PI, 3, _nf);
  init_delta_const(_nf);
  init_Wilson(mH);
  init_Wilson_mt(mH);
  //
  //muM = mH;
  muFpdf = _muF_mu_ratio * mH;
  //muF = _muF_mu_ratio * mH;
  //muR = _muR_muF_ratio * muF;
  mut = _muR_muF_ratio * _muF_mu_ratio * mH;
  //initResConsts(mH);
}





  double FO_1_log_PL(double z, double D(int,double), double muF_mH_ratio) {
    return C1log_pointlike_gg(z, D, _nf, muF_mH_ratio);
  }
  double FO_2_log_PL(double z, double D(int,double), double muF_mH_ratio) {
    return C2log_pointlike_gg(z, D, _nf, muF_mH_ratio);
  }
  double FO_3_log_PL(double z, double D(int,double), double muF_mH_ratio) {
    return C3log_pointlike_gg(z, D, _nf, muF_mH_ratio);
  }









// Fixed-order (NLO)
///////////////////////////////////////
double FO_1_log(double z, double D(int,double)) {
  //double lf = -2.*log(_muF_mu_ratio);
  //return FO_1_log_PL(z,D) + 2*CA*D(0,z)*lf;
  return C1log_pointlike_gg(z,D,_nf,_muF_mu_ratio);
}
double FO_1_delta_const() {
  double W1 = (finite_mt ? _W1t : _W1);
  double c = _cd1 + W1;
  //double lrf = 2.*log(_muR_muF_ratio);
  //c += 2*b0*lrf;
  return c;
}
// finite mt R_gg term
double Rgg_zspace(double z, double mu) {
  double z0 = 1e-6; // switches to asymptotic behaviour for z<z0
  double c10 = sxc->c10(1);
  double rgg = (z<z0 ? (2*c10+2.*log(z))/z : Rgg(z,_mH2_,qmasses,_integration_v_));
  return rgg;
}
double Rqg_zspace(double z, double mu) {
  double z0 = 1e-6; // switches to asymptotic behaviour for z<z0
  double c10 = sxc->c10(1);
  double rqg = (z<z0 ? (c10 + log(z))/z : Rqg(z,_mH2_,qmasses,_integration_v_));
  return CF*rqg;
}
double Rqqbar_zspace(double z, double mu) {
  // This has no numerical instabilities (no numerical integration)
  return Rqqbar(z,_mH2_,qmasses);
}
double FO_1_nonlog(double z, double mu) {
  double lf = -2.*log(_muF_mu_ratio);
  double smallx = (finite_mt ? CA*Rgg_zspace(z, mu) : - 11./2.*pow(1.-z,3.)/z );
  double nonlog = CA * 2.*( (1./z-2.+z*(1.-z))*(2.*log(1.-z)-log(z)) -log(z)/(1.-z) ) + smallx;
  nonlog += 2.*CA*(1./z-2.+z-z*z)*lf;
  return nonlog;
}
double FO_1_qg(double z, double mu) {
  double lf = -2.*log(_muF_mu_ratio);
  double Pgq = CF*(1.+pow(1.-z,2))/z/2.;
  double smallx = (finite_mt ? Rqg_zspace(z, mu) : 2./3.*z - pow(1.-z,2)/z );
  double res = Pgq*(2.*log(1-z)-log(z)+lf) + smallx;
  return res;
}
double FO_1_qqb(double z, double mu) {
  double res = (finite_mt ? Rqqbar_zspace(z, mu) : 32./27.*pow(1-z,3)/z );
  return res;
}
//




// Fixed-order (NNLO)
///////////////////////////////////////
double FO_2_log(double z, double D(int,double)) {
  double W1 = (finite_mt ? _W1t : _W1);
  double Clog = C2log_pointlike_gg(z,D,_nf,_muF_mu_ratio) + W1*C1log_pointlike_gg(z,D,_nf,_muF_mu_ratio);
  return Clog;
}
double FO_2_delta_const(double mH) {
  //double lf = -2.*log(_muF_mu_ratio);
  double W1 = (finite_mt ? _W1t : _W1);
  double W2 = (finite_mt ? _W2t : _W2);
  double c = _cd2 + W2 + _cd1*W1;
  //c += -18.*ZETA2*lf*lf + (-3./4. + 171./2.*ZETA3 + 1./4.*_nf - (33.-2.*_nf)/2.*ZETA2 - W1*(33.-2.*_nf)/12.)*lf;
  //
  //if(finite_mt) c = gg_nnlo_mt_delta(mu, mtop, _muF_mu_ratio); // no longer needed (checked)
  //
  //double lrf = 2.*log(_muR_muF_ratio);
  //double muRF = _muR_muF_ratio;  _muR_muF_ratio = 1.;
  //c += (2*b1+3*b0*FO_1_delta_const())*lrf + 3*b0*b0*lrf*lrf;
  //_muR_muF_ratio = muRF;
  return c;
}
double FO_2_smallx(double z, string channel, int mode, int truncation, int nlogs = 2) {
  //if(mode == -1 && z>0.5) return 0;
  double c10 = sxc->c10(_muF_mu_ratio);
  double c20 = sxc->c20(_muF_mu_ratio);
  double c11 = sxc->c11(_muF_mu_ratio);
  //cout << c10 << "   " << c20 << "   " << c11 << endl;
  if(channel=="gg") {}
  else if(channel=="qg") {
    c10 /= 2;
    c20 /= 2;
  }
  else if(channel=="qq") {
    c10 = 0;
    c20 = 0;
  } else {
    cout << "Error - wrong small-x channel" << endl;
    exit(98);
  }
  double c1 = 2*c10;
  double c2 = 2.*c20+c11;
  //
  double le00 = (-11.*CA+2.*_nf*(2.*CF/CA-1.))/12.;
  double le11 = _nf*(26.*CF-23.*CA)/36.;
  //
  double sub1=0, sub0=0;
  if(mode == 1) {
    for(int i=0; i<=truncation; i++) {
      sub1 += log1z_softexp[i] * pow(z-1,i);
      sub0 += pow(1-z,i);
    }
  }
  double hard = 0;
  if(nlogs>0) hard += c2*e01*e01 * ( -log(z)/z - sub1 );
  if(nlogs>1) hard += (2*c2*e01*le00 + c1*le11 - 2*c20*e01*b0) * ( 1./z - sub0 );
  //
  if(mode <= 0) hard *= pow(1-z,truncation+1);
  //if(mode ==-1) hard *= pow(1-2*z,truncation+1);
  if(channel=="qg") hard *= CF/CA;
  if(channel=="qq") hard *= pow(CF/CA,2);
  //
  return hard;
}
// sigma(z) = z C(z)
// C(z) = [D(z)]_+ + nonlog(z)
// sigma(z) = [D(z)]_+ - (1-z)D(z) + z*nonlog(z)
//          = [D(z)]_+ + reg(z)
// nonlog(z) = [ reg(z) + (1-z)D(z) ] /z
double FO_2_nonlog(double z, double mH) {
  double nonlog;
  if(!finite_mt) nonlog = C2reg_pointlike_gg(z, _nf, _muF_mu_ratio) + _W1 * C1reg_pointlike_gg(z, _nf, _muF_mu_ratio);
  else {
    if(mH<2*mtop) {
      nonlog = gg_nnlo_mt_reg(z, mH, mtop, _muF_mu_ratio) + FO_2_smallx(z,"gg",SXmode,Kact,2-SXlogs);
    } else {
      nonlog = gg_nnlo_infmt_reg(z, mH, mtop, _muF_mu_ratio);
      //nonlog = ( gg_nnlo_infmt_reg(z, mH, mtop, _muF_mu_ratio) + (1.-z)*gg_nnlo_infmt_distr(z, mH, mtop, _muF_mu_ratio) ) /z;
      //finite_mt = false;  nonlog = ( C2reg_pointlike(z, _nf, _muF_mu_ratio) + (1.-z)*FO_2_log(z,Dk) ) /z;  finite_mt = true;
      double nonlog1 = CA * 2.*( (1./z-2.+z*(1.-z))*(2.*log(1.-z)-log(z)) -log(z)/(1.-z) );
      nonlog += (_W1t-_W1)*nonlog1;
    }
  }
  //double lrf = 2.*log(_muR_muF_ratio);
  //nonlog += 3*b0*FO_1_nonlog(z,mH)*lrf;
  if(isPseudoscalar) {
    nonlog += deltaPseudoC2reg_gg(z,_nf);
  }
  return nonlog;
}
//
double FO_2_qg(double z, double mu) {
  double res;
  if(!finite_mt) res = C2reg_pointlike_qg(z, _nf, _muF_mu_ratio) + _W1*C1reg_pointlike_qg(z, _nf, _muF_mu_ratio);
  else {
    if(mu<2*mtop)
      res = qg_nnlo_mt(z, mu, mtop, _muF_mu_ratio) + FO_2_smallx(z,"qg",SXmode,Kact,2-SXlogs);
    else {
      res = qg_nnlo_infmt(z, mu, mtop, _muF_mu_ratio);
      double Pgq = CF*(1.+pow(1.-z,2))/z/2.;
      double res1 = Pgq*(2.*log(1-z)-log(z));
      res += (_W1t-_W1)*res1;
    }
  }
  if(isPseudoscalar) {
    res += deltaPseudoC2reg_qg(z,_nf);
  }
  return res;
}
double FO_2_qqb(double z, double mu) {
  double res;
  if(!finite_mt) res = C2reg_pointlike_qqb(z, _nf, _muF_mu_ratio) + _W1*C1reg_pointlike_qqb(z, _nf, _muF_mu_ratio) ;
  else {
    if(mu<2*mtop)
      res = qqb_nnlo_mt(z, mu, mtop, _muF_mu_ratio) + FO_2_smallx(z,"qq",SXmode,Kact,2-SXlogs);
    else {
      res = qqb_nnlo_infmt(z, mu, mtop, _muF_mu_ratio);
    }
  }
  return res;
  if(isPseudoscalar) {
    res += deltaPseudoC2reg_qqb(z,_nf);
  }
}
double FO_2_qq(double z, double mu) {
  double res;
  if(!finite_mt) res = C2reg_pointlike_qq(z, _nf, _muF_mu_ratio);
  else {
    if(mu<2*mtop)
      res = qq_nnlo_mt(z, mu, mtop, _muF_mu_ratio) + FO_2_smallx(z,"qq",SXmode,Kact,2-SXlogs);
    else {
      res = qq_nnlo_infmt(z, mu, mtop, _muF_mu_ratio);
    }
  }
  if(isPseudoscalar) {
    res += deltaPseudoC2reg_qq(z,_nf);
  }
  return res;
}
double FO_2_qqp(double z, double mu) {
  double res;
  if(!finite_mt) res = C2reg_pointlike_qqp(z, _nf, _muF_mu_ratio);
  else {
    if(mu<2*mtop)
      res = qqp_nnlo_mt(z, mu, mtop, _muF_mu_ratio) + FO_2_smallx(z,"qq",SXmode,Kact,2-SXlogs);
    else {
      res = qqp_nnlo_infmt(z, mu, mtop, _muF_mu_ratio);
    }
  }
  if(isPseudoscalar) {
    res += deltaPseudoC2reg_qqp(z,_nf);
  }
  return res;
}





  

// Fixed-order (NNNLO)
///////////////////////////////////////
double FO_3_log(double z, double D(int,double)) {
  double W1 = (finite_mt ? _W1t : _W1 );
  double W2 = (finite_mt ? _W2t : _W2 );
  double C3log = C3log_pointlike_gg(z,D,_nf,_muF_mu_ratio) + W1*C2log_pointlike_gg(z,D,_nf,_muF_mu_ratio) + W2*C1log_pointlike_gg(z,D,_nf,_muF_mu_ratio);
  return C3log;
}
double FO_3_delta_const(double mH) {
  //double lf = -2.*log(_muF_mu_ratio);
  // not yet known in the full theory => use pointlike value.
  //double W1 = (finite_mt ? _W1t : _W1 );
  //double W2 = (finite_mt ? _W2t : _W2 );
  //double W3 = (finite_mt ? _W3t : _W3 );
  double W1 = _W1;
  double W2 = _W2;
  double W3 = _W3;
  //
  double c = W3 + W2*_cd1 + W1*_cd2 + _cd3;
  return c;
}
double FO_3_smallx(double z, string channel, int mode, int truncation, int nlogs = 3) {
  if(mode == -1 && z>0.5) return 0;
  double c10 = sxc->c10(_muF_mu_ratio);
  double c20 = sxc->c20(_muF_mu_ratio);
  double c11 = sxc->c11(_muF_mu_ratio);
  double c30 = sxc->c30(_muF_mu_ratio);
  double c21 = sxc->c21(_muF_mu_ratio);
  if(channel=="gg") {}
  else if(channel=="qg") {
    c10 /= 2;
    c20 /= 2;
    c30 /= 2;
  }
  else if(channel=="qq") {
    c10 = 0;
    c20 = 0;
    c30 = 0;
  } else {
    cout << "Error - wrong small-x channel" << endl;
    exit(98);
  }
  if(MSbar) c30 += 8./3.*ZETA3;
  double c1 = 2*c10;
  double c2 = 2.*c20+c11;
  double c3 = 2.*(c30+c21);
  //
  double le00 = (-11.*CA+2.*_nf*(2.*CF/CA-1.))/12.;
  double le11 = _nf*(26.*CF-23.*CA)/36.;
  double le22 = -CA*( (18.*ZETA2-71.)*(CA-2*CF)*_nf + CA*CA*(54*ZETA3+99.*ZETA2-395) ) /108.;
  double le10 = 19.6798969734649 + _nf*(4./9.*ZETA2-68./81.) + _nf*_nf*13./2187.;
  double le21 = 222.097 - 5.73888*_nf + 0.113797*_nf*_nf;
  //
  double sub2=0, sub1=0, sub0=0;
  if(mode == 1) {
    for(int i=0; i<=truncation; i++) {
      sub2 += log2z_softexp[i] * pow(z-1,i);
      sub1 += log1z_softexp[i] * pow(z-1,i);
      sub0 += pow(1-z,i);
    }
  }
  double hard = 0;
  if(nlogs>0) hard +=  c3*e01*e01*e01 * ( pow(log(z),2)/z - sub2 )/2;
  if(nlogs>1) hard += (c3*e01*e01*le00*3 -(3*c30+c21)*2*b0*e01*e01  +c1*le22 +2*c2*e01*le11)* ( -log(z)/z - sub1 );
  if(nlogs>2) hard += (c3*e01*le00*le00*3-(3*c30+c21)*4*b0*e01*le00 +c1*le21 +2*c2*(le00*le11+e01*le10) -4*c20*b0*le11 +4*c30*b0*b0*e01) *( 1./z - sub0 );
  //
  if(mode == 0) hard *= pow(1-z,truncation+1);
  if(mode ==-1) hard *= pow(1-2*z,truncation+1);
  if(channel=="qg") hard *= CF/CA;
  if(channel=="qq") hard *= pow(CF/CA,2);
  //
  return hard;
}
  int TruncationOrder = 37;
  int C3mode = 2;
  void SetN3LOmode(int mode, int truncationOrder) {
    C3mode = mode;
    TruncationOrder = truncationOrder;
  }
const bool PseudoN3LOA = false;
double FO_3_nonlog(double z, double mu) {
  // there is no factorization of virtual contributions beyond the soft limit, so there is no point in using the "finit-mt" versions of the Wilson coeffs
  //double W1 = (finite_mt ? _W1t : _W1 );
  //double W2 = (finite_mt ? _W2t : _W2 );
  double W1 = _W1;
  double W2 = _W2;
  double res = 0;
  if(finite_mt) {
    res = W1*C2reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + W2*C1reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    res += FO_3_smallx(z, "gg", SXmode, TruncationOrder, 3-SXlogs);
    if(SXmode==-1) {
      res = W1*C2reg_pointlike_gg_softmidexp(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_gg_softmidexp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_gg_softmidexp(z, _nf, _muF_mu_ratio);
      res += FO_3_smallx(z, "gg", SXmode, 50, 3-SXlogs);
    }
  } else {
    //if(C3mode==0) res = W1*C2reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + W2*C1reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    //if(C3mode==1) res = W1*(C2reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder)+C2reg_pointlike_gg_smallx(z, _nf, TruncationOrder)) + W2*(C1reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder)+C1reg_pointlike_gg_smallx(z, _nf, TruncationOrder)) + C3reg_pointlike_gg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    //if(C3mode>=2)
    res = W1*C2reg_pointlike_gg(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_gg(z, _nf, _muF_mu_ratio) + C3reg_pointlike_gg(z, _nf, _muF_mu_ratio);
    if(PseudoN3LOA && isPseudoscalar) res += 0.5*deltaPseudoC2reg_gg(z,_nf);
  }
  return res;
}
double FO_3_qg(double z, double mu) {
  //double W1 = (finite_mt ? _W1t : _W1 );
  //double W2 = (finite_mt ? _W2t : _W2 );
  double W1 = _W1;
  double W2 = _W2;
  double res = 0;
  if(finite_mt) {
    res = W1*C2reg_pointlike_qg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + W2*C1reg_pointlike_qg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_qg_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    res += FO_3_smallx(z, "qg", SXmode, TruncationOrder, 3-SXlogs);
    if(SXmode==-1) {
      res = W1*C2reg_pointlike_qg_softmidexp(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_qg_softmidexp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qg_softmidexp(z, _nf, _muF_mu_ratio);
      res += FO_3_smallx(z, "qg", SXmode, 50, 3-SXlogs);
    }
  } else {
    res = W1*C2reg_pointlike_qg(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_qg(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qg(z, _nf, _muF_mu_ratio);
    if(PseudoN3LOA && isPseudoscalar) res += 0.5*deltaPseudoC2reg_qg(z,_nf);
  }
  return res;
}
double FO_3_qqb(double z, double mu) {
  //double W1 = (finite_mt ? _W1t : _W1 );
  //double W2 = (finite_mt ? _W2t : _W2 );
  double W1 = _W1;
  double W2 = _W2;
  double res = 0;
  if(finite_mt) {
    res = W1*C2reg_pointlike_qqb_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + W2*C1reg_pointlike_qqb_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_qqb_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    res += FO_3_smallx(z, "qq", SXmode, TruncationOrder, 3-SXlogs);
    if(SXmode==-1) {
      res = W1*C2reg_pointlike_qqb_softmidexp(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_qqb_softmidexp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qqb_softmidexp(z, _nf, _muF_mu_ratio);
      res += FO_3_smallx(z, "qq", SXmode, 50, 3-SXlogs);
    }
  } else {
    res = W1*C2reg_pointlike_qqb(z, _nf, _muF_mu_ratio) + W2*C1reg_pointlike_qqb(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qqb(z, _nf, _muF_mu_ratio);
    if(PseudoN3LOA && isPseudoscalar) res += 0.5*deltaPseudoC2reg_qqb(z,_nf);
  }
  return res;
}
double FO_3_qq(double z, double mu) {
  //double W1 = (finite_mt ? _W1t : _W1 );
  double W1 = _W1;
  double res = 0;
  if(finite_mt) {
    res = W1*C2reg_pointlike_qq_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_qq_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    res += FO_3_smallx(z, "qq", SXmode, TruncationOrder, 3-SXlogs);
    if(SXmode==-1) {
      res = W1*C2reg_pointlike_qq_softmidexp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qq_softmidexp(z, _nf, _muF_mu_ratio);
      res += FO_3_smallx(z, "qq", SXmode, 50, 3-SXlogs);
    }
  } else {
    res = W1*C2reg_pointlike_qq(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qq(z, _nf, _muF_mu_ratio);
    if(PseudoN3LOA && isPseudoscalar) res += 0.5*deltaPseudoC2reg_qq(z,_nf);
  }
  return res;
}
double FO_3_qqp(double z, double mu) {
  //double W1 = (finite_mt ? _W1t : _W1 );
  double W1 = _W1;
  double res = 0;
  if(finite_mt) {
    res = W1*C2reg_pointlike_qqp_softexp(z, _nf, _muF_mu_ratio, TruncationOrder) + C3reg_pointlike_qqp_softexp(z, _nf, _muF_mu_ratio, TruncationOrder);
    res += FO_3_smallx(z, "qq", SXmode, TruncationOrder, 3-SXlogs);
    if(SXmode==-1) {
      res = W1*C2reg_pointlike_qqp_softmidexp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qqp_softmidexp(z, _nf, _muF_mu_ratio);
      res += FO_3_smallx(z, "qq", SXmode, 50, 3-SXlogs);
    }
  } else {
    res = W1*C2reg_pointlike_qqp(z, _nf, _muF_mu_ratio) + C3reg_pointlike_qqp(z, _nf, _muF_mu_ratio);
    if(PseudoN3LOA && isPseudoscalar) res += 0.5*deltaPseudoC2reg_qqp(z,_nf);
  }
  return res;
}











double FO_log(double z, double D(int,double)) {
  double res = 0;
  if(order == 1) res = FO_1_log(z,D);
  if(order == 2) res = FO_2_log(z,D);
  if(order == 3) res = FO_3_log(z,D);
  return res;
}
double FO_distr_gg(double z, double mu) {
  double res = 0;
  res = FO_log(z,Dk);
  return res;
}
double FO_reg_gg(double z, double mu) {
  double res = 0;
  if(smallxRes && computeSmallxRes) {
    res =  SxRES_gg(_alphas, z, _muF_mu_ratio, order);
  } else {
    if(order == 1) res = FO_1_nonlog(z,mu);
    if(order == 2) res = FO_2_nonlog(z,mu);
    if(order == 3) res = FO_3_nonlog(z,mu);
  }
  return res;
}
double FO_delta_gg(double mu) {
  double res = 0;
  if(order == 1) res = FO_1_delta_const();
  if(order == 2) res = FO_2_delta_const(mu);
  if(order == 3) res = FO_3_delta_const(mu);
  return res;
}
double FO_qg(double z, double mu) {
  double res = 0;
  if(smallxRes && computeSmallxRes) {
    res =  SxRES_qg(_alphas, z, _muF_mu_ratio, order);
  } else {
    if(order == 1) res = FO_1_qg(z, mu);
    if(order == 2) res = FO_2_qg(z, mu);
    if(order == 3) res = FO_3_qg(z, mu);
  }
  return res;
}
double FO_qqb(double z, double mu) {
  double res = 0;
  if(smallxRes && computeSmallxRes) {
    res = SxRES_qq(_alphas, z, _muF_mu_ratio, order);
  } else {
    if(order == 1) res = FO_1_qqb(z, mu);
    if(order == 2) res = FO_2_qqb(z, mu);
    if(order == 3) res = FO_3_qqb(z, mu);
  }
  return res;
}
double FO_qq(double z, double mu) {
  double res = 0;
  if(smallxRes && computeSmallxRes) {
    res = SxRES_qq(_alphas, z, _muF_mu_ratio, order);
  } else {
    if(order == 2) res = FO_2_qq(z, mu);
    if(order == 3) res = FO_3_qq(z, mu);
  }
  return res;
}
double FO_qqp(double z, double mu) {
  double res = 0;
  if(smallxRes && computeSmallxRes) {
    res = SxRES_qq(_alphas, z, _muF_mu_ratio, order);
  } else {
    if(order == 2) res = FO_2_qqp(z, mu);
    if(order == 3) res = FO_3_qqp(z, mu);
  }
  return res;
}











// gg channel
// integrand of regular part
double gg_reg_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = lumi(lumi_gg, t2, tau/z, muFpdf) * FO_reg_gg(z, mu) /z;
  return res * (1-tau);
}
// integrand of distributional part
double gg_distr_int(double t1, double t2, double tau, double mu) {
  if(computeSmallxRes) return 0;
  double lumtau = lumi(lumi_gg, t2, tau, muFpdf);
  double res = 0;//lumtau * FO_delta_gg(mu);  // delta term
  if(order>0) {
    double z = 1.-(1.-tau)*t1;
    res += (1-tau) * ( lumi(lumi_gg, t2, tau/z, muFpdf)/z - lumtau ) * FO_distr_gg(z, mu);
    res += - lumtau * FO_distr_gg(t1*tau, mu) * tau; // add residual integral form 0 to tau
  }
  return res;
}
//
// non-gg channels contributions
double qg_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = lumi(lumi_qg, t2, tau/z, muFpdf) * FO_qg(z, mu) /z;
  return res * (1-tau);
}
//
  /*
double qqb_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = lumi(lumi_qqb, t2, tau/z, muFpdf) * FO_qqb(z, mu) /z;
  return res * (1-tau);
}
//
double qq_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = lumi(lumi_qq, t2, tau/z, muFpdf) * FO_qq(z, mu) /z;
  return res * (1-tau);
}
//
double qqp_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = lumi(lumi_qqp, t2, tau/z, muFpdf) * FO_qqp(z, mu) /z;
  return res * (1-tau);
}
  */
double quark_channels_int(double t1, double t2, double tau, double mu) {
  double z = 1.-(1.-tau)*t1;
  double res = 0;
  res += lumi(lumi_qqb, t2, tau/z, muFpdf) * FO_qqb(z, mu) /z;
  if(order>1) {
    res += lumi(lumi_qq , t2, tau/z, muFpdf) * FO_qq (z, mu) /z;
    res += lumi(lumi_qqp, t2, tau/z, muFpdf) * FO_qqp(z, mu) /z;
  }
  return res * (1-tau);
  // version with non-singlet basis
  double qqb = FO_qqb(z, mu);
  double qq  = FO_qq (z, mu);
  double qqp = FO_qqp(z, mu);
  double FO_SS = (qq + qqb + 2*(_nf-1.)*qqp)/(2.*_nf);
  double FO_TT = (qq + qqb - 2*qqp)/2.;
  double FO_VV = (qq - qqb)/2.;
  res += lumi(lumi_SS, t2, tau/z, muFpdf) * FO_SS /z;
  res += lumi(lumi_VV, t2, tau/z, muFpdf) * FO_VV /z;
  res += lumi(lumi_TT, t2, tau/z, muFpdf) * FO_TT /z;
  return res * (1-tau);
}




struct int1dp { double tau, mu; };
bool __all_channels__ = false;
int integrandCUBA(const int *ndim, const double x[], const int *ncomp, double res[], void *pars) {
  int1dp p = *(int1dp *) pars;
  double tau = p.tau;
  double mu  = p.mu;
  if(finite_mt && order==1 && !computeSmallxRes) _integration_v_ = x[2];
  res[0] = 0;
  if(_gg_channel_) {
    res[0] += gg_reg_int(x[0],x[1],tau,mu);
    res[0] += gg_distr_int(x[0],x[1],tau,mu);
  }
  if(__all_channels__) {
    if(_qg_channel_) res[0] += qg_int(x[0],x[1],tau,mu);
    if(_qq_channel_) res[0] += quark_channels_int(x[0],x[1],tau,mu);
    //res[0] += qqb_int(x[0],x[1],tau,mu);
    //if(order > 1) res[0] += qq_int(x[0],x[1],tau,mu) + qqp_int(x[0],x[1],tau,mu);
  }
  return 0;
}
double integral2d(double tau, double mu, double *err = NULL) {
  double res = 0;
  int1dp p = { tau, mu };
  //
  double cres[1], error[1], prob[1];
  double epsrel=INTEGRAL_PRECISION, epsabs=1.e-10;
  int last = 4;
  int verbose = 0;
  int nregions, neval, fail;
  //system("export CUBACORES=8");
  int dim = 2;
  if(finite_mt && order==1 && !computeSmallxRes) dim = 3;
  if(finite_mt && computeSmallxRes && order<=2) epsrel *= 5;
  if(finite_mt && computeSmallxRes && order==3) epsrel *= 20;
  Cuhre(dim, 1, integrandCUBA, &p, 1,
	epsrel, epsabs, verbose | last,
	0, 500000, 9,
	NULL, NULL,
	&nregions, &neval, &fail, cres, error, prob);
  // Suave(2, 1, integrandCUBA, &p,
  // 	epsrel, epsabs, verbose | last, 0,
  // 	0, 500000, 1000, 25.,
  // 	&nregions, &neval, &fail, cres, error, prob);
  // Vegas(2, 1, integrandCUBA, &p,
  // 	epsrel, epsabs, verbose | last, 0,
  // 	0, 500000, 1000, 500, 1000,
  // 	0, NULL,
  // 	&neval, &fail, cres, error, prob);
  //  
  if(fail == -1) {
    cout << "\033[0;31m  ERROR: dimension out of range.\033[0m" << endl;
  } else if(fail>0) {
    cout << "\033[0;31m  ERROR: precision not reached.\033[0m" << endl;
  }
  double lum_tau = Lum(tau, muFpdf);
  res += cres[0]/lum_tau;
  if(err != NULL) *err = error[0]/lum_tau;
  return res;
}


double xs(double mH, double sqrts, double *err = NULL) {
  double tau = mH*mH/sqrts/sqrts;
  double res = integral2d(tau, mH, err);
  if(_gg_channel_ && !computeSmallxRes) res += FO_delta_gg(mH);
  return res;
}



  /*
dcomplex brn_q_ep0(dcomplex xin) { 
  dcomplex fx, Ax;
  double x = real(xin);
  if(x>=1) fx = pow(asin(1/sqrt(x)),2);
  else fx = -1./4.*pow(log((1.+sqrt(1.-x))/(1.-sqrt(1.-x)))-M_PI*I,2);
  //Ax = 3.*x/2.*(1.+(1.-x)*fx);
  Ax = x*fx;
  return Ax;
}
dcomplex sborn_mass(double mH, double mt, double mb, double mc){
  dcomplex res = 0;
  for(int i=0; i<number_of_quarks; i++) {
    dcomplex Wq = 4.*m[i]*m[i]/mH/mH;
    res += brn_q_ep0(Wq);
  }
  return res/M_PI;
}
  */
//dcomplex brn_q_ep0(dcomplex x) {
//  dcomplex logx = log(x);
//  dcomplex HPL2 = 0.5*logx*logx;
//  return -16.*(1.+x)*(1.+x)/x*HPL2 + 32.*(x-1.)*(x-1.)/x;
//}
dcomplex brn_q_ep0(dcomplex Wq) {
  dcomplex res = 0;
  if(isPseudoscalar) {
    dcomplex fx;
    double x = real(Wq);
    if(x>=1) fx = pow(asin(1/sqrt(x)),2);
    else fx = -1./4.*pow(log((1.+sqrt(1.-x))/(1.-sqrt(1.-x)))-M_PI*I,2);
    res = x*fx*3./2.;
  } else {
    dcomplex x = -Wq/pow(1.+sqrt(1.-Wq),2);
    dcomplex logx = log(x);
    dcomplex HPL2 = 0.5*logx*logx;
    res = -16.*(1.+x)*(1.+x)/x*HPL2 + 32.*(x-1.)*(x-1.)/x;
    res *= Wq*Wq /64.*(-3./4.);
  }
  return res;
}
dcomplex sborn_mass(dcomplex mH2, vector<double> m) {
  dcomplex res = 0;
  for(int i=0; i<m.size(); i++) {
    dcomplex Wq = 4.*m[i]*m[i]/mH2;
    //dcomplex Xq = -Wq/pow(1.+sqrt(1.-Wq),2);
    //res += brn_q_ep0(Xq)*Wq*Wq /64.*(-3./4.);
    res += brn_q_ep0(Wq);
  }
  return res;
}
double xs_LO_prefactor_EFT() {
  double pref = 35.0309/M_PI/M_PI; // = Gf/pi/sqrt(2)/288 with the Gf in pb
  return pref;
}
double EFTrescaling(double mH, double mt, double mb, double mc) {
  vector<double> m;
  if(mt>0) m.push_back(mt);
  if(mb>0) m.push_back(mb);
  if(mc>0) m.push_back(mc);
  dcomplex mH2 = mH*mH + I*1e-10; // analytic continuation
  return norm(sborn_mass(mH2,m));
}
double xs_LO_prefactor(double mH) {
  return xs_LO_prefactor_EFT() * EFTrescaling(mH,mtop,mbot,mcha);
}
  double EFTrescaling_deriv(double mH, double mt, double mb, double mc, int id, double &second_deriv, double &third_deriv) {
    double m0[3] = {mt, mb, mc};
    double m [3] = {mt, mb, mc};
    double eps = 0.2; // cross-checked with Mathematica for the top - agreement within 6 digits even for the third derivative
    const int np = 9;
    double res[np];
    for(int i=0; i<np; i++) {
      m [id] = m0[id] + eps*(i-(np-1)/2);
      res[i] = EFTrescaling(mH,m[0],m[1],m[2]);
    }
    //double fd_coeff[np] = { -1./60., 3./20., -3./4., 0, 3./4., -3./20., 1./60. };
    //double sd_coeff[np] = { 1./90., -3./20., 3./2., -49./18., 3./2., -3./20., 1./90. };
    //double td_coeff[np] = { 1./8., -1., 13./8., 0, -13./8., 1., -1./8. };
    double fd_coeff[np] = { 1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280. };
    double sd_coeff[np] = { -1./560., 8./315., -1./5., 8./5., -205./72., 8./5., -1./5., 8./315., -1./560. };
    double td_coeff[np] = { -7./240., 3./10., -169./120., 61./30., 0, -61./30., 169./120., -3./10., 7./240. };
    double first_deriv = 0;
    second_deriv = 0;
    third_deriv = 0;
    for(int i=0; i<np; i++) {
      first_deriv  += fd_coeff[i] * res[i];
      second_deriv += sd_coeff[i] * res[i];
      third_deriv  += td_coeff[i] * res[i];
    }
    first_deriv  /= eps;
    second_deriv /= eps*eps;
    third_deriv  /= eps*eps*eps;
    return first_deriv;
  }

  double EFT_rescaling_MSbar_NLO(double mH) {
    double fd, sd, td;
    fd = EFTrescaling_deriv(mH, mtop, mbot, mcha, 0, sd, td);
    return fd * mtop * alphas(mut)/M_PI * (-4./3.-2*log(mut/mtmt));
  }






  void init_compute(double mH) {
    if(isPseudoscalar && finite_mt) {
      cout << "Error: pseudoscalar at finite top mass not available." << endl;
      exit(-37);
    }
    if(finite_mt) init_nnlo_mt_coeff();
    qmasses.clear();
    if(mtop>0) qmasses.push_back(mtop);
    if(mbot>0) qmasses.push_back(mbot);
    if(mcha>0) qmasses.push_back(mcha);
    _mH2_ = mH*mH + I*1e-10; // analytic continuation
    _alphas = alphas(mH*_muF_mu_ratio*_muR_muF_ratio); // as(muR)
    init_consts(mH);
    return;
  }



bool firstcall = true;
void compute(double mH, double sqrts, double LO, double NLO, double NNLO, double *rexact) {
  init_compute(mH);
  double as = _alphas;
  //cout << "alpha_s = "<< as << endl;
  double api = as/M_PI;
  double fact = pow(api,order) * LO;
  //double tau = mH*mH/sqrts/sqrts;
  //
  if(!quiet && firstcall) cout << "sigma = sigma0 * [ 1 + as K1 + as^2 K2 + as^3 K3 + ...]" << endl;
  firstcall = false;
  //
  string channel = (__all_channels__ ? "all channels:" : "gg channel:");
  if(!quiet) cout << channel << endl;
  string sassig[4] = { "  sig0 = ",
		       "  as   K1 = ",
		       "  as^2 K2 = ",
		       "  as^3 K3 = " };
  double err;
  double lrf = 2.*log(_muR_muF_ratio);
  double mur_dep_term = 0; // using terms of previous orders *at* mur
  if(!_gg_channel_) LO = 0;
  NNLO -= NLO;
  NLO -= LO;
  if(order == 1) mur_dep_term = api*LO*2*b0*lrf;
  if(order == 2) mur_dep_term = api*(api*2*b1*LO + 3*b0*NLO) *lrf - api*api*3*b0*b0*LO*lrf*lrf;
  //if(order == 2) mur_dep_term = api*(api*2*b1*LO + 3*b0*NLO) *lr - api*api*3*b0*b0*LO*lr*lr;
  if(order == 3) mur_dep_term = api*(api*api*2*b2*LO + api*3*b1*NLO + 4*b0*NNLO) *lrf
		   - api*api*(api*7*b0*b1*LO + 6*b0*b0*NLO) *lrf*lrf + api*api*api*4*b0*b0*b0*LO*lrf*lrf*lrf;
  //cout << order << "  " << mur_dep_term << endl;
  //
  *rexact = xs(mH, sqrts, &err) * fact + mur_dep_term;
  err *= fact;
  if(!quiet) cout << sassig[order] << *rexact/LO << " +/- " << err/LO
		  << "  (exact)" << endl;
}

  double computeSxRes(double mH, double sqrts, double LO, int ord, double inputalphas) {
#ifndef withHELL
    return 0;
#endif
    //if(order==3) return 0;
    order = ord;
    init_compute(mH);
    if(inputalphas > 0) _alphas = inputalphas;
    double err;
    computeSmallxRes = true;
    double res = xs(mH, sqrts, &err) * LO;
    computeSmallxRes = false;
    return res;
  }



  double GetZmass() {
    return mZ;
  }


};

double ChebLum(double z, int LumChannel) {
  return ggHiggs::Lumij(z, ggHiggs::muFpdf, LumChannel);
}


