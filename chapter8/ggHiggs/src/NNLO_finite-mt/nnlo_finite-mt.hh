#pragma once
#include "../math/complex.hh"
#include "nnlo_pars.hh"

namespace ggHiggs {

// z-space
void init_nnlo_mt_coeff();
double gg_nnlo_mt_delta(double mH, double mt, double muF_mH_ratio=1.);
double gg_nnlo_mt_distr(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nnlo_mt_reg(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nnlo_infmt_distr(double z, double mH, double mt, double muF_mH_ratio=1.);
double gg_nnlo_infmt_reg(double z, double mH, double mt, double muF_mH_ratio=1.);
int    gg_nnlo_mt_delta_ser(double mH, double mt, double *res, double muF_mH_ratio=1.);
// other channels
double qg_nnlo_mt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qq_nnlo_mt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qqp_nnlo_mt(double z, double mH, double mt, double muF_mH_ratio=1.);
double qqb_nnlo_mt(double z, double mH, double mt, double muF_mH_ratio=1.);
double qg_nnlo_infmt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qq_nnlo_infmt (double z, double mH, double mt, double muF_mH_ratio=1.);
double qqp_nnlo_infmt(double z, double mH, double mt, double muF_mH_ratio=1.);
double qqb_nnlo_infmt(double z, double mH, double mt, double muF_mH_ratio=1.);


// N-space
dcomplex gg_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex gg_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
// other channels
dcomplex qg_nnlo_mt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qq_nnlo_mt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqp_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqb_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qg_nnlo_infmt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qq_nnlo_infmt_N (dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqp_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);
dcomplex qqb_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio=1.);


// mathcing at small z (or N)
double   smallx_match_nnlo_z(double z, double c2);
dcomplex smallx_match_nnlo(dcomplex N, double c2);
  //dcomplex smallx_match_der(dcomplex N, double c2);

double smallx_gg_match_nnlo_z(double z, double c2, double mHmt);

};
