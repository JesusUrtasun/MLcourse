dcomplex gg_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res =  -2.1363609499084208e+02*psi(1.0000000000000000e+00,N)-nl*lmH*z2+-4.5000000000000001e-01*(z2*z2)+6.6666666666666663e-01*lmtop*nl+2.3750000000000000e+00*lmtop+5.4000000000000000e+01*psi(N)*psi(2.0000000000000000e+00,N)+8.5433157316424515e-13*(xht*xht)*( -2.2004269056000000e+11*psi(1.0000000000000000e+00,N)+-1.0425515520000000e+10*lmtop*nl+9.1197388800000000e+10*lmtop+8.8369711242750000e+12*z3+2.2004269056000000e+11*pow(psi(N),2.0)+-1.3268828160000000e+10*nl*z2+-3.2256000000000000e+05*( 9.7486740000000000e+06*z3+-4.5827100000000000e+05*(lmH*lmH)+2.4996600000000000e+06*lmH*z2+2.7774000000000000e+04*nl*(lmH*lmH)+-1.1109600000000000e+05*nl*z2+8.6408000000000000e+04*nl+-5.0158620000000000e+06*lmH+1.8330840000000000e+06*z2+9.2580000000000000e+04*nl*lmH-1.8701160000000000e+06)*( psi(N)+5.7721566490153287e-01)+1.8701160448000000e+10*nl+5.5157760000000000e+07*( 6.3180000000000000e+03*z3+-2.9700000000000000e+02*(lmH*lmH)+1.6200000000000000e+03*lmH*z2+1.8000000000000000e+01*nl*(lmH*lmH)+-7.2000000000000000e+01*nl*z2+5.6000000000000000e+01*nl+-2.3940000000000000e+03*lmH+1.1880000000000000e+03*z2+6.0000000000000000e+01*nl*lmH-1.2120000000000000e+03)*( psi(N)+5.7721566490153287e-01)+2.5402417587662531e+11*psi(N)+6.3221760000000000e+07*( 4.4226000000000000e+04*z3+-2.0790000000000000e+03*(lmH*lmH)+1.1340000000000000e+04*lmH*z2+1.2600000000000000e+02*nl*(lmH*lmH)+-5.0400000000000000e+02*nl*z2+3.9200000000000000e+02*nl+-2.0022000000000000e+04*lmH+8.3160000000000000e+03*z2+4.2000000000000000e+02*nl*lmH-8.4840000000000000e+03)*( psi(N)+5.7721566490153287e-01)+9.6555110400000000e+09*lmH+3.1805016117687683e+11*z2+4.3132185600000000e+09*nl*lmH-9.5511826068370684e+12)+-4.1250000000000000e+01*z3+2.1363609499084208e+02*pow(psi(N),2.0)+-8.3119055745820731e+01*psi(1.0000000000000000e+00,N)*psi(N)+-4.1559527872910365e+01*( psi(1.0000000000000000e+00,N)-pow(psi(N),2.0))*psi(N)+-8.0375514403292187e-08*xht*( 1.8800640000000000e+07*psi(1.0000000000000000e+00,N)+1.2139200000000000e+06*lmtop*nl+-1.0172160000000000e+07*lmtop+-4.2956572500000000e+08*z3+-1.8800640000000000e+07*pow(psi(N),2.0)+9.6768000000000000e+05*nl*z2+-3.7281280000000000e+06*nl+-4.0320000000000000e+04*( 6.3180000000000000e+03*z3+-2.9700000000000000e+02*(lmH*lmH)+1.6200000000000000e+03*lmH*z2+1.8000000000000000e+01*nl*(lmH*lmH)+-7.2000000000000000e+01*nl*z2+5.6000000000000000e+01*nl+-2.3940000000000000e+03*lmH+1.1880000000000000e+03*z2+6.0000000000000000e+01*nl*lmH-1.2120000000000000e+03)*( psi(N)+5.7721566490153287e-01)+-2.1704047836348709e+07*psi(N)+5.7600000000000000e+03*( 4.4226000000000000e+04*z3+-2.0790000000000000e+03*(lmH*lmH)+1.1340000000000000e+04*lmH*z2+1.2600000000000000e+02*nl*(lmH*lmH)+-5.0400000000000000e+02*nl*z2+3.9200000000000000e+02*nl+-2.0022000000000000e+04*lmH+8.3160000000000000e+03*z2+4.2000000000000000e+02*nl*lmH-8.4840000000000000e+03)*( psi(N)+5.7721566490153287e-01)+1.5552000000000000e+06*lmH+-2.5948209327368494e+07*z2+-6.9168000000000000e+05*nl*lmH+4.3718221858399421e+08)+-2.7777777777777776e-02*( 6.0000000000000000e+00*psi(1.0000000000000000e+00,N)+-6.0000000000000000e+00*pow(psi(N),2.0)+-6.9265879788183939e+00*psi(N)-1.1868671943935670e+01)*( 1.0800000000000000e+02*(lmH*lmH)+-1.0000000000000000e+01*nl+9.9000000000000000e+01*lmH+-2.7000000000000000e+02*z2+-6.0000000000000000e+00*nl*lmH+3.9900000000000000e+02)+1.6500000000000000e+01*lmH*z2+8.3333333333333337e-01*z3*nl+-8.5500000000000000e+01*z3*lmH+-1.6666666666666667e+00*nl*z2+-1.8000000000000000e+01*(lmH*lmH)*z2+-8.2569444444444446e+00*nl+-2.7777777777777776e-02*( 6.3180000000000000e+03*z3+-2.9700000000000000e+02*(lmH*lmH)+1.6200000000000000e+03*lmH*z2+1.8000000000000000e+01*nl*(lmH*lmH)+-7.2000000000000000e+01*nl*z2+5.6000000000000000e+01*nl+-2.3940000000000000e+03*lmH+1.1880000000000000e+03*z2+6.0000000000000000e+01*nl*lmH-1.2120000000000000e+03)*( psi(N)+5.7721566490153287e-01)+-1.8000000000000000e+01*psi(3.0000000000000000e+00,N)+5.5555555555555552e-02*( 2.0000000000000000e+00*nl+-1.0800000000000000e+02*lmH-3.3000000000000000e+01)*( 1.0389881968227591e+01*psi(1.0000000000000000e+00,N)+-1.0389881968227591e+01*pow(psi(N),2.0)+1.2000000000000000e+01*psi(1.0000000000000000e+00,N)*psi(N)+6.0000000000000000e+00*( psi(1.0000000000000000e+00,N)-pow(psi(N),2.0))*psi(N)+-3.5606015831807014e+01*psi(N)+-6.0000000000000000e+00*psi(2.0000000000000000e+00,N)-3.2669246738911909e+01)+3.9203096086694285e+02*psi(N)+1.3500000000000000e+01*lmH+5.4000000000000000e+01*psi(1.0000000000000000e+00,N)*( psi(1.0000000000000000e+00,N)-pow(psi(N),2.0))+4.1559527872910365e+01*psi(2.0000000000000000e+00,N)+-1.8000000000000000e+01*( 2.0000000000000000e+00*psi(1.0000000000000000e+00,N)*psi(N)+( psi(1.0000000000000000e+00,N)-pow(psi(N),2.0))*psi(N)-psi(2.0000000000000000e+00,N))*psi(N)+6.6500000000000000e+01*z2+-1.8333333333333333e+00*nl*lmH+5.0326625573468311e+02;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_gg__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,0);
    res += coeff_nnlo_gg__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_gg__1_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__1_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,1);
    res += coeff_nnlo_gg__1_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__1_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__1_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__1_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__1_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_gg__2_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_gg__2_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,2);
    res += coeff_nnlo_gg__2_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_gg__2_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_gg__2_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_gg__2_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_gg__2_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
  }
  return res;
}

dcomplex qg_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qg__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,0);
    res += coeff_nnlo_qg__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qg__1_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__1_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,1);
    res += coeff_nnlo_qg__1_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__1_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__1_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__1_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__1_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qg__2_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__2_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,2);
    res += coeff_nnlo_qg__2_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__2_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__2_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__2_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__2_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qg__3_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qg__3_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,3);
    res += coeff_nnlo_qg__3_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qg__3_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qg__3_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qg__3_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qg__3_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
  }
  return res;
}

dcomplex qq_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qq__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,0);
    res += coeff_nnlo_qq__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qq__1_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__1_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,1);
    res += coeff_nnlo_qq__1_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__1_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__1_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__1_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__1_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qq__2_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__2_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,2);
    res += coeff_nnlo_qq__2_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__2_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__2_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__2_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__2_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qq__3_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qq__3_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,3);
    res += coeff_nnlo_qq__3_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qq__3_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qq__3_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qq__3_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qq__3_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
  }
  return res;
}

dcomplex qqp_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qqp_0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,0);
    res += coeff_nnlo_qqp_0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqp_1_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_1_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,1);
    res += coeff_nnlo_qqp_1_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_1_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_1_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_1_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_1_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqp_2_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_2_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,2);
    res += coeff_nnlo_qqp_2_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_2_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_2_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_2_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_2_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqp_3_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqp_3_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,3);
    res += coeff_nnlo_qqp_3_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqp_3_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqp_3_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqp_3_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqp_3_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
  }
  return res;
}

dcomplex qqb_nnlo_mt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double xht = mH*mH/mt/mt;
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qqb_0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,0);
    res += coeff_nnlo_qqb_0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,0);
    res += coeff_nnlo_qqb_1_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_1_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,1);
    res += coeff_nnlo_qqb_1_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_1_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_1_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_1_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_1_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,1);
    res += coeff_nnlo_qqb_2_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_2_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,2);
    res += coeff_nnlo_qqb_2_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_2_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_2_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_2_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_2_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,2);
    res += coeff_nnlo_qqb_3_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqb_3_0_0_1[k] * mell_L(k,0,N) * pow(lmH,0) * pow(lmtop,1) * pow(xht,3);
    res += coeff_nnlo_qqb_3_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqb_3_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqb_3_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqb_3_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0) * pow(lmtop,0) * pow(xht,3);
    res += coeff_nnlo_qqb_3_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1) * pow(lmtop,0) * pow(xht,3);
  }
  return res;
}

dcomplex gg_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double lmH = 2.*log(muF_mH_ratio);
  double lmtop = 2.*log(muF_mH_ratio*mH/mt);
  dcomplex res =  12.0*( 12.0*zeta(3.0)+2.0*(Pi*Pi)*Euler+( 6.0*(Euler*Euler)+(Pi*Pi))*Euler)*psi(N)-nl*lmH*z2+-(9.0/20.0)*(z2*z2)+(2.0/3.0)*lmtop*nl+(19.0/8.0)*lmtop+54.0*psi(N)*psi(2.0,N)+(6.0/5.0)*((Pi*Pi)*(Pi*Pi))+(3.0/2.0)*( 6.0*(Euler*Euler)+(Pi*Pi))*(Pi*Pi)+-(165.0/4.0)*z3+3.0*( 12.0*zeta(3.0)+2.0*(Pi*Pi)*Euler+( 6.0*(Euler*Euler)+(Pi*Pi))*Euler)*Euler+(33.0/2.0)*lmH*z2+-18.0*( 6.0*(Euler*Euler)+(Pi*Pi))*( psi(1.0,N)-pow(psi(N),2.0))+(5.0/6.0)*z3*nl-( 108.0*(lmH*lmH)+-10.0*nl+99.0*lmH+-270.0*z2+-6.0*nl*lmH+399.0)*( 6.0*psi(1.0,N)+-6.0*(Euler*Euler)+-12.0*psi(N)*Euler+-6.0*pow(psi(N),2.0)-(Pi*Pi))/36.0+-(171.0/2.0)*z3*lmH+-(5.0/3.0)*nl*z2+-18.0*(lmH*lmH)*z2+-(1189.0/144.0)*nl-( psi(N)+Euler)*( 6318.0*z3+-297.0*(lmH*lmH)+1620.0*lmH*z2+18.0*nl*(lmH*lmH)+-72.0*nl*z2+56.0*nl+-2394.0*lmH+1188.0*z2+60.0*nl*lmH-1212.0)/36.0+-18.0*psi(3.0,N)+-72.0*( 2.0*psi(1.0,N)*psi(N)+( psi(1.0,N)-pow(psi(N),2.0))*psi(N)-psi(2.0,N))*Euler+( 18.0*( psi(1.0,N)-pow(psi(N),2.0))*Euler+12.0*psi(1.0,N)*psi(N)+6.0*( psi(1.0,N)-pow(psi(N),2.0))*psi(N)+-12.0*zeta(3.0)+-2.0*(Pi*Pi)*Euler-( 6.0*(Euler*Euler)+(Pi*Pi))*Euler+-6.0*psi(2.0,N)+-3.0*( 6.0*(Euler*Euler)+(Pi*Pi))*psi(N))*( 2.0*nl+-108.0*lmH-33.0)/18.0+(27.0/2.0)*lmH+54.0*psi(1.0,N)*( psi(1.0,N)-pow(psi(N),2.0))+-18.0*( 2.0*psi(1.0,N)*psi(N)+( psi(1.0,N)-pow(psi(N),2.0))*psi(N)-psi(2.0,N))*psi(N)+(133.0/2.0)*z2+108.0*zeta(3.0)*Euler+-(11.0/6.0)*nl*lmH+(11399.0/144.0);
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_gg__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0);
    res += coeff_nnlo_gg__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1);
    res += coeff_nnlo_gg__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2);
    res += coeff_nnlo_gg__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0);
    res += coeff_nnlo_gg__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1);
    res += coeff_nnlo_gg__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2);
    res += coeff_nnlo_gg__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0);
    res += coeff_nnlo_gg__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1);
    res += coeff_nnlo_gg__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0);
  }
  return res;
}

dcomplex qg_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double lmH = 2.*log(muF_mH_ratio);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qg__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0);
    res += coeff_nnlo_qg__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1);
    res += coeff_nnlo_qg__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2);
    res += coeff_nnlo_qg__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0);
    res += coeff_nnlo_qg__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1);
    res += coeff_nnlo_qg__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2);
    res += coeff_nnlo_qg__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0);
    res += coeff_nnlo_qg__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1);
    res += coeff_nnlo_qg__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0);
  }
  return res;
}

dcomplex qq_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double lmH = 2.*log(muF_mH_ratio);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qq__0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0);
    res += coeff_nnlo_qq__0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1);
    res += coeff_nnlo_qq__0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2);
    res += coeff_nnlo_qq__0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0);
    res += coeff_nnlo_qq__0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1);
    res += coeff_nnlo_qq__0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2);
    res += coeff_nnlo_qq__0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0);
    res += coeff_nnlo_qq__0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1);
    res += coeff_nnlo_qq__0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0);
  }
  return res;
}

dcomplex qqp_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double lmH = 2.*log(muF_mH_ratio);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qqp_0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0);
    res += coeff_nnlo_qqp_0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1);
    res += coeff_nnlo_qqp_0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2);
    res += coeff_nnlo_qqp_0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0);
    res += coeff_nnlo_qqp_0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1);
    res += coeff_nnlo_qqp_0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2);
    res += coeff_nnlo_qqp_0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0);
    res += coeff_nnlo_qqp_0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1);
    res += coeff_nnlo_qqp_0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0);
  }
  return res;
}

dcomplex qqb_nnlo_infmt_N(dcomplex N, double mH, double mt, double muF_mH_ratio) {
  double lmH = 2.*log(muF_mH_ratio);
  dcomplex res = 0;
  for(int k=0; k<14; k++) {
    res += coeff_nnlo_qqb_0_0_0_0[k] * mell_L(k,0,N) * pow(lmH,0);
    res += coeff_nnlo_qqb_0_0_1_0[k] * mell_L(k,0,N) * pow(lmH,1);
    res += coeff_nnlo_qqb_0_0_2_0[k] * mell_L(k,0,N) * pow(lmH,2);
    res += coeff_nnlo_qqb_0_1_0_0[k] * mell_L(k,1,N) * pow(lmH,0);
    res += coeff_nnlo_qqb_0_1_1_0[k] * mell_L(k,1,N) * pow(lmH,1);
    res += coeff_nnlo_qqb_0_1_2_0[k] * mell_L(k,1,N) * pow(lmH,2);
    res += coeff_nnlo_qqb_0_2_0_0[k] * mell_L(k,2,N) * pow(lmH,0);
    res += coeff_nnlo_qqb_0_2_1_0[k] * mell_L(k,2,N) * pow(lmH,1);
    res += coeff_nnlo_qqb_0_3_0_0[k] * mell_L(k,3,N) * pow(lmH,0);
  }
  return res;
}

