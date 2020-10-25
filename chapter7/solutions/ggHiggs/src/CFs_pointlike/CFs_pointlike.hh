#pragma once
#include "../parameters.hh"
#include "../math/complex.hh"
#include <iostream>
#include <cstdlib>

// Writing C_ij = (1 + as W1 + as^2 W2 + as^3 W3) * (delta_ig delta_jg + as Cp1_ij + as^2 Cp2_ij + as^3 Cp3_ij)
// these are the Cpn_ij terms


// order as
extern double C1reg_pointlike_gg (double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C1reg_pointlike_qg (double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C1reg_pointlike_qqb(double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
//
extern double C1reg_pointlike_gg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C1reg_pointlike_qg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C1reg_pointlike_qqb_softexp(double z, int nf, double muF_mH_ratio, int truncation=37);
//
extern double C1reg_pointlike_gg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C1reg_pointlike_qg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C1reg_pointlike_qqb_midexp(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C1reg_pointlike_gg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C1reg_pointlike_qg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C1reg_pointlike_qqb_midexpplain(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C1reg_pointlike_gg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C1reg_pointlike_qg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C1reg_pointlike_qqb_softmidexp(double z, int nf, double muF_mH_ratio);
//
extern double C1reg_pointlike_gg_smallx  (double z, int nf, int truncation=37);



// order as^2
extern double C2reg_pointlike_gg (double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C2reg_pointlike_qg (double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C2reg_pointlike_qqb(double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C2reg_pointlike_qq (double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
extern double C2reg_pointlike_qqp(double z, int nf, double muf_mu_ratio=1., double mur_muf_ratio=1.);
//
extern double C2reg_pointlike_gg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C2reg_pointlike_qg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C2reg_pointlike_qqb_softexp(double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C2reg_pointlike_qq_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C2reg_pointlike_qqp_softexp(double z, int nf, double muF_mH_ratio, int truncation=37);
//
extern double C2reg_pointlike_gg_smallx (double z, int nf, int truncation=37);
extern double C2reg_pointlike_qg_smallx (double z, int nf, int truncation=37);
extern double C2reg_pointlike_qqb_smallx(double z, int nf, int truncation=37);
extern double C2reg_pointlike_qq_smallx (double z, int nf, int truncation=37);
extern double C2reg_pointlike_qqp_smallx(double z, int nf, int truncation=37);
//
extern double C2reg_pointlike_gg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qqb_midexp(double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qq_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qqp_midexp(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C2reg_pointlike_gg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qqb_midexpplain(double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qq_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C2reg_pointlike_qqp_midexpplain(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C2reg_pointlike_gg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C2reg_pointlike_qg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C2reg_pointlike_qqb_softmidexp(double z, int nf, double muF_mH_ratio);
extern double C2reg_pointlike_qq_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C2reg_pointlike_qqp_softmidexp(double z, int nf, double muF_mH_ratio);
//
extern dcomplex C2reg_pointlike_gg_softexp (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qg_softexp (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qqb_softexp(dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qq_softexp (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qqp_softexp(dcomplex N, int nf, int truncation=37);
//
extern dcomplex C2reg_pointlike_gg_smallx (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qg_smallx (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qqb_smallx(dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qq_smallx (dcomplex N, int nf, int truncation=37);
extern dcomplex C2reg_pointlike_qqp_smallx(dcomplex N, int nf, int truncation=37);
//
extern dcomplex C2_pointlike_gg (dcomplex N);
extern dcomplex C2_pointlike_qg (dcomplex N);
extern dcomplex C2_pointlike_qqb(dcomplex N);
extern dcomplex C2_pointlike_qq (dcomplex N);
extern dcomplex C2_pointlike_qqp(dcomplex N);
// Pseudo scalar
extern double deltaPseudoC2reg_gg (double z, int nf);
extern double deltaPseudoC2reg_qg (double z, int nf);
extern double deltaPseudoC2reg_qqb(double z, int nf);
extern double deltaPseudoC2reg_qq (double z, int nf);
extern double deltaPseudoC2reg_qqp(double z, int nf);



// order as^3
extern double C3reg_pointlike_gg (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qg (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qqb(double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qq (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qqp(double z, int nf, double muF_mH_ratio);
//
extern double C3reg_pointlike_gg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37, double power=0);
extern double C3reg_pointlike_qg_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C3reg_pointlike_qqb_softexp(double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C3reg_pointlike_qq_softexp (double z, int nf, double muF_mH_ratio, int truncation=37);
extern double C3reg_pointlike_qqp_softexp(double z, int nf, double muF_mH_ratio, int truncation=37);
//
extern double C3reg_pointlike_gg_heexp (double z, int nf, double muF_mH_ratio, int truncation=38);
extern double C3reg_pointlike_qg_heexp (double z, int nf, double muF_mH_ratio, int truncation=38);
extern double C3reg_pointlike_qqb_heexp(double z, int nf, double muF_mH_ratio, int truncation=38);
extern double C3reg_pointlike_qq_heexp (double z, int nf, double muF_mH_ratio, int truncation=38);
extern double C3reg_pointlike_qqp_heexp(double z, int nf, double muF_mH_ratio, int truncation=38);
//
extern double C3reg_pointlike_gg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qg_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qqb_midexp(double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qq_midexp (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qqp_midexp(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C3reg_pointlike_gg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qg_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qqb_midexpplain(double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qq_midexpplain (double z, int nf, double muF_mH_ratio, int truncation=50);
extern double C3reg_pointlike_qqp_midexpplain(double z, int nf, double muF_mH_ratio, int truncation=50);
//
extern double C3reg_pointlike_gg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qg_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qqb_softmidexp(double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qq_softmidexp (double z, int nf, double muF_mH_ratio);
extern double C3reg_pointlike_qqp_softmidexp(double z, int nf, double muF_mH_ratio);
//
extern double C3reg_pointlike_gg_smallx(double z, int truncation=37);



template <class T>
T C1log_pointlike_gg(T z, T D(int,T), int nf, double muf_mu_ratio=1.) {
  double lf = -2.*log(muf_mu_ratio);
  return 4.*CA*D(1,z) + 2*CA*D(0,z)*lf;
}
template <class T>
T C2log_pointlike_gg(T z, T D(int,T), int nf, double muf_mu_ratio=1.) {
  double lf = -2.*log(muf_mu_ratio);
  return ( 72. * D(3,z)
	   + (2.*nf - 33.) * D(2,z)
	   + (67. - 90.*ZETA2 - 10./3.*nf) * D(1,z)
	   + (-101./3. + 33.*ZETA2 + 351./2.*ZETA3 + nf*(14./9.-2.*ZETA2)) * D(0,z)
	   )
    +    ( 108*lf * D(2,z)
	   + (36.*lf*lf - (33.-2.*nf)*lf) * D(1,z)
	   + (-(33.-2.*nf)/4.*lf*lf + (67./2.-5./3.*nf-5.*CA*CA*ZETA2)*lf) * D(0,z) );
}
template <class T>
T C3log_pointlike_gg(T z, T D(int,T), int nf, double muf_mu_ratio=1.) {
  double lf = -2.*log(muf_mu_ratio);
  //if(lf!=0) { std::cout << "Error! not yet implemented!" << std::endl; std::exit(9);}
  double CA3 = CA*CA*CA;
  return ( 8*CA3 * D(5,z)
	   + (-110./9.*CA3 + 20/9.*CA*CA*nf) * D(4,z)
	   + (CA3*(925./27. - 56*ZETA2) - 164./27.*CA*CA*nf + 4./27.*CA*nf*nf) * D(3,z)
	   + (CA3 * (-1051./27. + 187./3*ZETA2 + 181*ZETA3)
	      + CA*CA*nf * (457./2./27. - 34./3.*ZETA2)
	      + 0.5*CA*CF*nf - 10./27.*CA*nf*nf ) * D(2,z)
	   + (CA3 * (30569./8./81. - 152./3.*ZETA2 - 352./3.*ZETA3 - 154./5.*ZETA2*ZETA2)
	      + CA*CA*nf * (-4211./4./81. + 94./9.*ZETA2 + 46./3.*ZETA3)
	      - CA*CF*nf * (63./8.-6.*ZETA3)
	      + CA*nf*nf * (25./81.-4./9.*ZETA2) ) * D(1,z)
	   + (CA3 * (-297029./32./729. + 8563./4./81.*ZETA2 + 8941./4./27.*ZETA3 + 253./4./15.*ZETA2*ZETA2
		     - 725./2./3.*ZETA2*ZETA3 + 186.*ZETA5)
	      + CA*CA*nf * (31313./16./729. - 2173./4./81.*ZETA2 - 475./4./9.*ZETA3 - 17./2./15*ZETA2*ZETA2)
	      + CA*CF*nf * (1711./32./27. - 0.5*ZETA2 - 19./2./9.*ZETA3 - 1./5.*ZETA2*ZETA2)
	      - CA*nf*nf * (58./729. - 10./27.*ZETA2 - 5./27.*ZETA3) ) * D(0,z)
	   )
    // scale dep terms copied from ihixs 2.0 - they assume that the Wilson coefficient is computed at muF
    +    ( 540.*lf * D(4,z)
	   +( 432*lf*lf - (660-40*nf)*lf ) * D(3,z)
	   +( 108*lf*lf*lf + (-891./2.+27.*nf)*lf*lf + (2775./2. -378.*6*ZETA2 +2./3.*nf*nf -82.*nf)*lf ) * D(2,z)
	   +( (6*nf-99)*lf*lf*lf + (1971./4. +1./3.*nf*nf -162*6*ZETA2 -31*nf)*lf*lf
	      - (10./9.*nf*nf + (-130./3.+13.*6*ZETA2)*nf + 60*(51./8.-19./24.*nf) -429*3*ZETA2 -4860*ZETA3 +1257./2.)*lf ) * D(1,z)
	   +( (1./18.*nf*nf -11./6.*nf -108.*ZETA2 + 121./8.)*lf*lf*lf
	      + (-5./18.*nf*nf + (13./6.-21./2.*ZETA2)*nf + 945.*ZETA3 + 693./4.*ZETA2 - 161./8. - 27.*(51./8.-19./24.*nf))*lf*lf
	      - ((2./3.*ZETA2 -25./54.)*nf*nf + (5345./72. -81.*ZETA3 -47.*ZETA2)*nf - 30569./48. + 114.*6*ZETA2 + 1584.*ZETA3 + 2079./5.*ZETA2*ZETA2)*lf ) * D(0,z)
	   );
}


