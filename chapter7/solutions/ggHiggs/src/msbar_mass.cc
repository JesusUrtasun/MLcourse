#include "msbar_mass.hh"
#include <cmath>
#include "parameters.hh"
#include "ggHiggs.hh"


double beta0(int nf) { return (33.-2.*nf)/12.; }
double beta1(int nf) { return (153.-19.*nf)/24.; }
double beta2(int nf) { return (2857. - 5033./9.*nf + 325./27.*nf*nf)/128.; }
double beta3(int nf) { return (149753./6. + (1093.*nf*nf*nf)/729. + 3564*ZETA3
			       + nf*nf*(50065./162. + (6472.*ZETA3)/81.)
			       - nf*(1078361./162. + (6508.*ZETA3)/27.))/256.; }



// routine to compute ms-bar bottom mass
double runmass(double mass0, double api0, double apif, unsigned int order, int nf) {
  /*
   evaluates the running of the MS-bar quark mass
   by expanding the equation
   
   m(mu) = m(mu0) * exp( \int_a0^af dx gammam(x)/x/beta(x) )
   
   in terms of alpha_s. The results agree with RunDec.m.
   
   Input:
   ------
   mass0  :  m(mu0)
   api0   :  alpha_s(mu0)/pi
   apif   :  alpha_s(muf)/pi
   nf     :  number of flavors
   nloop  :  order of calculation (nloop=1..4)

   Output:
   -------
   massout:  m(muf)
  */
  //double nloop = order+1; // old way
  double nloop = order;
  //if (nloop == 0) return mass0;
  //
  double bb0 = beta0(nf);
  double bb1 = beta1(nf)/bb0;
  double bb2 = beta2(nf)/bb0;
  double bb3 = beta3(nf)/bb0;
  //
  double gamma0 = 1.;
  double gamma1 = (202./3. - (20.*nf)/9.)/16.;
  double gamma2 = (1249. - (140*nf*nf)/81.
		   + 2*nf*(-556./27.- 48.*ZETA3) + (8*nf*(-46. + 48.*ZETA3))/9.)/64.;
  double gamma3 = (28413.91975308642 + (135680*ZETA3)/27.
		   + nf*nf*nf*(-1.3662551440329218 + (64*ZETA3)/27.)
		   + nf*nf*(21.57201646090535 - (16.*pow(M_PI,4))/27. + (800*ZETA3)/9.)
		   - 8800*ZETA5 + nf*(-3397.1481481481483 + (88*pow(M_PI,4))/9.
				       - (34192*ZETA3)/9. + (18400*ZETA5)/9.))/256.;
  //
  double cc0 = gamma0/bb0;
  double cc1 = gamma1/bb0;
  double cc2 = gamma2/bb0;
  double cc3 = gamma3/bb0;
  //
  double cfunc1 = 1.;
  double cfunc2 = cc1 - bb1*cc0;
  double cfunc3 = 1./2.*(pow(cc1-bb1*cc0,2) + cc2 - bb1*cc1 + bb1*bb1*cc0 - bb2*cc0);
  double cfunc4 = (1./6.*pow(cc1 - bb1*cc0,3) + 1./2.*(cc1 - bb1*cc0)*(cc2 - bb1*cc1 + bb1*bb1*cc0 - bb2*cc0)
		   + 1./3.*(cc3 - bb1*cc2 + bb1*bb1*cc1 - bb2*cc1 - bb1*bb1*bb1*cc0 + 2*bb1*bb2*cc0 - bb3*cc0));
  //
  if (nloop < 4) {
    cfunc4 = 0.;
    if (nloop < 3) {
      cfunc3 = 0.;
      if (nloop < 2) {
	cfunc2 = 0.;
	if (nloop < 1) 
	  cfunc1 = 0.;
      }
    }
  }
  //
  double cfuncmu0 = cfunc1 + cfunc2*api0 + cfunc3*api0*api0 + cfunc4*api0*api0*api0;
  double cfuncmuf = cfunc1 + cfunc2*apif + cfunc3*apif*apif + cfunc4*apif*apif*apif;
  //cout << cfuncmu0 << endl;
  //cout << cfuncmuf << endl;
  //
  return mass0*pow(apif/api0,cc0)*cfuncmuf/cfuncmu0;
}


double pole_to_MSbar(double mbPole, double alphas, int order) {
  double res = mbPole;
  double api = alphas/M_PI;
  double nf = 5;
  if(order == 1) res = mbPole*(1-4./3.*api);
  if(order == 2) res = mbPole*(1-4./3.*api + api*api*(ZETA3/6.-ZETA2/3.*(2.*log(2.)+7)-2393./288.+nf/3.*(ZETA2+71./48.)));
  return res;
}

double MSbar_to_pole(double mbmb, double alphas, int order) { // converts MSbar to pole mass, using value of mb(mb)
  double res = mbmb;
  double api = alphas/M_PI;
  double nf = 5;
  if(order == 1) res = mbmb*(1+4./3.*api);
  if(order == 2) res = mbmb*(1+4./3.*api + api*api*(307./32.-(71.*nf)/144.+2*ZETA2 +(2*log(2.)*ZETA2)/3.-(nf*ZETA2)/3.-ZETA3/6.) );
  return res;
}
