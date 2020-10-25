dcomplex gg_nlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  dcomplex res =  6.0*pow(psi(N),2.0)+6.0*(Euler*Euler)+12.0*psi(N)*Euler+-6.0*psi(1.0,N)+2.0*(Pi*Pi)+_W1t;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_gg_0_0[k] * mell_L(k,0,N) * pow(xht,0);
    res += coeff_nlo_gg_0_1[k] * mell_L(k,1,N) * pow(xht,0);
    res += coeff_nlo_gg_1_0[k] * mell_L(k,0,N) * pow(xht,1);
    res += coeff_nlo_gg_1_1[k] * mell_L(k,1,N) * pow(xht,1);
    res += coeff_nlo_gg_2_0[k] * mell_L(k,0,N) * pow(xht,2);
    res += coeff_nlo_gg_2_1[k] * mell_L(k,1,N) * pow(xht,2);
    res += coeff_nlo_gg_3_0[k] * mell_L(k,0,N) * pow(xht,3);
    res += coeff_nlo_gg_3_1[k] * mell_L(k,1,N) * pow(xht,3);
  }
  return res;
}

dcomplex qg_nlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  dcomplex res = 0;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_qg_0_0[k] * mell_L(k,0,N) * pow(xht,0);
    res += coeff_nlo_qg_0_1[k] * mell_L(k,1,N) * pow(xht,0);
    res += coeff_nlo_qg_1_0[k] * mell_L(k,0,N) * pow(xht,1);
    res += coeff_nlo_qg_1_1[k] * mell_L(k,1,N) * pow(xht,1);
    res += coeff_nlo_qg_2_0[k] * mell_L(k,0,N) * pow(xht,2);
    res += coeff_nlo_qg_2_1[k] * mell_L(k,1,N) * pow(xht,2);
    res += coeff_nlo_qg_3_0[k] * mell_L(k,0,N) * pow(xht,3);
    res += coeff_nlo_qg_3_1[k] * mell_L(k,1,N) * pow(xht,3);
  }
  return res;
}

dcomplex qqb_nlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  dcomplex res = 0;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_qq_0_0[k] * mell_L(k,0,N) * pow(xht,0);
    res += coeff_nlo_qq_0_1[k] * mell_L(k,1,N) * pow(xht,0);
    res += coeff_nlo_qq_1_0[k] * mell_L(k,0,N) * pow(xht,1);
    res += coeff_nlo_qq_1_1[k] * mell_L(k,1,N) * pow(xht,1);
    res += coeff_nlo_qq_2_0[k] * mell_L(k,0,N) * pow(xht,2);
    res += coeff_nlo_qq_2_1[k] * mell_L(k,1,N) * pow(xht,2);
    res += coeff_nlo_qq_3_0[k] * mell_L(k,0,N) * pow(xht,3);
    res += coeff_nlo_qq_3_1[k] * mell_L(k,1,N) * pow(xht,3);
  }
  return res;
}

dcomplex gg_nlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  dcomplex res =  6.0*pow(psi(N),2.0)+6.0*(Euler*Euler)+12.0*psi(N)*Euler+-6.0*psi(1.0,N)+2.0*(Pi*Pi)+_W1t;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_gg_0_0[k] * mell_L(k,0,N);
    res += coeff_nlo_gg_0_1[k] * mell_L(k,1,N);
  }
  return res;
}

dcomplex qg_nlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  dcomplex res = 0;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_qg_0_0[k] * mell_L(k,0,N);
    res += coeff_nlo_qg_0_1[k] * mell_L(k,1,N);
  }
  return res;
}

dcomplex qqb_nlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  dcomplex res = 0;
  for(int k=0; k<27; k++) {
    res += coeff_nlo_qq_0_0[k] * mell_L(k,0,N);
    res += coeff_nlo_qq_0_1[k] * mell_L(k,1,N);
  }
  return res;
}

