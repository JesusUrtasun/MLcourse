#pragma once
#include "../math/complex.hh"
#include "nlo_pars.hh"

namespace ggHiggs {

// z-space
double gg_nlo_mt_delta(double mH, double mt, double muF_mH_ratio=1.);
double gg_nlo_mt_distr(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nlo_mt_reg(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nlo_infmt_distr(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nlo_infmt_reg(double z, double mH, double mt, double muF_mH_ratio=1.);
// other channels
double qg_nlo_mt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qqb_nlo_mt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qg_nlo_infmt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qqb_nlo_infmt (double z, double mH, double mt, double muF_mH_ratio=1.);


// N-space
dcomplex gg_nlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
  //dcomplex gg_nlo_infmt_der_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex gg_nlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
  //dcomplex gg_nlo_mt_der_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
// other channels
dcomplex qg_nlo_mt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqb_nlo_mt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qg_nlo_infmt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqb_nlo_infmt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);


// mathcing at small z (or N)
double   smallx_match_nlo_z(double z, double c1);
dcomplex smallx_match_nlo_N(dcomplex N, double c1);
  //dcomplex smallx_match_nlo_N_der(dcomplex N, double c1);

};
