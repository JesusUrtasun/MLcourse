#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "nlo_finite-mt.hh"
#include "../math/mellin_Dk.hh"

using namespace std;

extern dcomplex psi(dcomplex x);
extern dcomplex psi(int n, dcomplex x);

namespace ggHiggs {

const double Pi = M_PI;
const double Euler = 0.577215664901532860606512090082402431042159335939;
//const double zeta_2 = 1.6449340668482264366;
//const double zeta_3 = 1.2020569031595942854;
/*
double zeta(double x) {
  double res;
  if(x==2) res = zeta_2;
  else if(x==3) res = zeta_3;
  else exit(2);
  return res;
}
*/

extern double _W1t;

// z space
#include "nlo_mt_zspace.cc"


// N space ( expansion in (1-z) )
#include "nlo_mt_Nspace.cc"



double smallx_match_nlo_z(double z, double c1) {
  double res = 1./z;
  // j<Kact  is consistent with the other program
  for(int j=1; j<KactNLO; j++) {
    res -= pow(1-z,j);
  }
  return res*c1*3.;// /Pi/Pi;
}
dcomplex smallx_match_nlo_N(dcomplex N, double c1) {
  dcomplex res = 1./(N-1);
  // j<Kact  is consistent with the other program
  for(int j=1; j<KactNLO; j++) {
    for(int p=0; p<=j; p++) {
      res -= pow(-1.,p)*gsl_sf_choose(j,p)/(N+p);
    }
  }
  return res*c1*3.;// /Pi/Pi;
}
// derivative
// dcomplex smallx_match_nlo_N_der(dcomplex N, double c1) {
//   dcomplex res = -1./(N-1)/(N-1);
//   // j<Kact  is consistent with the other program
//   for(int j=1; j<KactNLO; j++) {
//     for(int p=0; p<=j; p++) {
//       res += pow(-1.,p)*gsl_sf_choose(j,p)/(N+p)/(N+p);
//     }
//   }
//   return res*c1*3.;// /Pi/Pi;
// }


};
