#pragma once

#include <cmath>
#include <vector>
#include "../math/complex.hh"

/*
  computation of the coefficients h_{k1,k2} of the expansion of the impact factor

  (1)   h_Higgs(M1,M2) = sum_{k1,k2>=0} h_{k1,k2} M1^k1 M2^k2

  where M1 and M2 are the Mellin conjugate variables of the two gluon off-shellness.
*/

class HiggsSmallX {
public:
  HiggsSmallX() { init(125, 173); };
  HiggsSmallX(double mH, double mt, double mb=-1, double mc=-1, double muFrat=1) { init(mH,mt,mb,mc,muFrat); };
  ~HiggsSmallX() { return; };
  void   SetFactorizationScale(double muFrat) { _muFrat=muFrat; };
  double GetFactorizationScale() { return _muFrat; };
  void init(double mH, double mt, double mb=-1, double mc=-1, double muFrat=1);
  // this function returns a vector A such that h_{k1,k2} = A[k1+k2][k1], up to k1+k2 <= maxorder
  std::vector<std::vector<double> > hCoeffs(int maxorder=3);
  double h1();
  double factorH();
  double h_kp_km(int kp, int km, double *err = NULL);
private:
  std::vector<dcomplex> Yi;
  double prefH;
  double _muFrat;
};



/*
// functions below should be removed (not needed in the main file)
extern double integrand2d_ch_var(double x, double y, int kp, int km, std::vector<dcomplex> Y);
//extern double integrand1d_ch_var(double x, void *p);
extern double integrand2d(double t, double y, int kp, int km, std::vector<dcomplex> Y);
//extern double integrand1d(double t, int kp);
extern double integrand_h1(double x, void *params);

extern dcomplex cA1    (double t, double y, dcomplex yt);
extern dcomplex cA1dt  (double t, double y, dcomplex yt);
extern dcomplex cA1dy  (double t, double y, dcomplex yt);
extern dcomplex cA1dtdy(double t, double y, dcomplex yt);

extern dcomplex cA1y0  (double t, dcomplex yt);
extern dcomplex cA1y0dt(double t, dcomplex yt);

extern dcomplex cA1y1  (double t, dcomplex yt);
extern dcomplex cA1y1dt(double t, dcomplex yt);

extern dcomplex cA3    (double t, double y, dcomplex yt);

extern dcomplex C0    (double t, double y, dcomplex yt); // to delete
extern dcomplex C0dt  (double t, double y, dcomplex yt); // to delete
extern dcomplex C0dy  (double t, double y, dcomplex yt); // to delete
extern dcomplex C0dtdy(double t, double y, dcomplex yt); // to delete
extern dcomplex C0dtdt(double t, double y, dcomplex yt); // to delete
extern dcomplex C0dydy(double t, double y, dcomplex yt); // to delete
extern dcomplex gC0dtdt(double t, double y, dcomplex yt); // to delete
extern dcomplex B0(double rho, dcomplex cyt); // to delete


//extern dcomplex ddB0(double rho, dcomplex cyt);

*/



