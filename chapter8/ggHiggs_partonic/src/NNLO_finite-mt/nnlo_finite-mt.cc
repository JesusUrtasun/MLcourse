#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "nnlo_finite-mt.hh"
#include "../math/mellin_Dk.hh"
#include "../math/special_functions.hh"

using namespace std;

extern dcomplex psi(dcomplex x);
extern dcomplex psi(int n, dcomplex x);

namespace ggHiggs {

const double Pi = M_PI;
const double Euler = 0.577215664901532860606512090082402431042159335939;
const double zeta_2 = 1.6449340668482264366;
const double zeta_3 = 1.2020569031595942854;
const double nl = 5.;

double zeta(double x) {
  double res;
  if(x==2) res = zeta_2;
  else if(x==3) res = zeta_3;
  else exit(2);
  return res;
}
const double z2=zeta(2.);
const double z3=zeta(3.);








  const double SEpower = 0;   // 0 means expansion of C(z), 1 means expansion of z*C(z) (as in Harlander) [SEpower = -a, where a is the variable in the HELL3 paper]
  double lastSEpower = 2147483647;


// z space
#include "nnlo_mt_zspace.cc"






// N space ( expansion in (1-z) )
#include "nnlo_mt_Nspace.cc"







double smallx_match_nnlo_z(double z, double c2) {
  double res = -log(z);
  // j<Kact  is consistent with the other program
  //for(int j=1; j<Kact; j++) {
  //  res -= pow(1-z,j)/(j+0.);
  //}
  return res*c2*9./z *pow(1-z,2)*pow(1-sqrt(z),4);//* pow(1-z, 10*z/(1-z));
}
dcomplex smallx_match_nnlo(dcomplex N, double c2) {
  dcomplex res = 1./N/N;
  // j<Kact  is consistent with the other program
  //for(int j=1; j<Kact; j++) {
  //  for(int p=0; p<=j; p++) {
  //    res -= pow(-1.,p)*gsl_sf_choose(j,p)/(N+p)/(j+0.);
  //  }
  //}
  return res*c2*9.;// /Pi/Pi;
}
// derivative
// dcomplex smallx_match_nnlo_der(dcomplex N, double c2) {
//   dcomplex res = -2./N/N/N;
//   // j<Kact  is consistent with the other program
//   for(int j=1; j<Kact; j++) {
//     for(int p=0; p<=j; p++) {
//       res += pow(-1.,p)*gsl_sf_choose(j,p)/(N+p)/(N+p)/(j+0.);
//     }
//   }
//   return res*c2*9.;// /Pi/Pi;
// }

  /*
  double sx_gg_0  [] = {-281.592, 232.34, 63.1424, 30.419, 18.4763, 12.7782, 9.61335, 7.67202, 6.39417, 5.50714, 4.86517, 4.38459, 4.01456, 3.7228, 3.48799, 3.29559};
  double sx_gg_1  [] = {-1.27652, 2.47533, 1.29972, 0.862577, 0.621373, 0.463071, 0.348483, 0.260167, 0.18908, 0.130018, 0.0797504, 0.0361533, -0.0022352, -0.0364589, -0.0672857, -0.0952957};
  double sx_gg_2_1[] = {-0.0691071, 0.0958929, 0.133929, 0.171964, 0.21, 0.248036, 0.286071, 0.324107, 0.362143, 0.400179, 0.438214, 0.47625, 0.514286, 0.552321, 0.590357, 0.628393};
  double sx_gg_2_0[] = {-0.641418, 0.985466, 0.556801, 0.443316, 0.401032, 0.384062, 0.378103, 0.377459, 0.379507, 0.382913, 0.386954, 0.391221, 0.395474, 0.39957, 0.403427, 0.406996};
  double smallx_gg_match_nnlo_z(double z, double c2, double mHmt) {
    double xht = mHmt*mHmt;
    double lHt = log(mHmt);
    double res = -9*c2*log(z);
    // j<Kact  is consistent with the other program
    for(int j=0; j<Kact; j++) {
      res -= pow(1-z,j) *( 0+sx_gg_0  [j]
			   + sx_gg_1  [j] * xht
			   + sx_gg_2_0[j] * xht*xht
			   + sx_gg_2_1[j] * xht*xht * lHt );
    }
    return res * pow(1-z, 10*z/(1-z));
  }
  /*/
  double sx_gg_0  [] = {-281.592, -49.2513, 13.891, 44.31, 62.7863, 75.5645, 85.1779, 92.8499, 99.244, 104.751, 109.616, 114.001, 118.016, 121.738, 125.226, 128.522};
  double sx_gg_1  [] = {-1.27652, 1.19882, 2.49853, 3.36111, 3.98248, 4.44556, 4.79404, 5.0542, 5.24328, 5.3733, 5.45305, 5.48921, 5.48697, 5.45051, 5.38323, 5.28793};
  double sx_gg_2_1[] = {-0.0691071, 0.0267857, 0.160714, 0.332679, 0.542679, 0.790714, 1.07679, 1.40089, 1.76304, 2.16321, 2.60143, 3.07768, 3.59196, 4.14429, 4.73464, 5.36304};
  double sx_gg_2_0[] = {-0.641418, 0.344047, 0.900849, 1.34416, 1.7452, 2.12926, 2.50736, 2.88482, 3.26433, 3.64724, 4.0342, 4.42542, 4.82089, 5.22046, 5.62389, 6.03088};
  //double sx_gg_0  [] = {0., 232.34, 295.483, 325.902, 344.378, 357.156, 366.769, 374.441, 380.836, 386.343, 391.208, 395.593, 399.607, 403.33, 406.818, 410.113};
  //double sx_gg_1  [] = {0.66731, 3.14264, 4.44236, 5.30494, 5.92631, 6.38938, 6.73787, 6.99803, 7.18711, 7.31713, 7.39688, 7.43303, 7.4308, 7.39434, 7.32705, 7.23176};
  //double sx_gg_2_1[] = {0.0578571, 0.15375, 0.287679, 0.459643, 0.669643, 0.917679, 1.20375, 1.52786, 1.89, 2.29018, 2.72839, 3.20464, 3.71893, 4.27125, 4.86161, 5.49};
  //double sx_gg_2_0[] = {-0.222933, 0.762532, 1.31933, 1.76265, 2.16368, 2.54774, 2.92585, 3.30331, 3.68281, 4.06573, 4.45268, 4.8439, 5.23937, 5.63895, 6.04237, 6.44937};
  double smallx_gg_match_nnlo_z(double z, double c2, double mHmt) {
    double xht = mHmt*mHmt;
    double lHt = log(mHmt);
    double res = 0;
    // j<Kact  is consistent with the other program
    double gg_0[Kact], gg_1[Kact], gg_2_1[Kact], gg_2_0[Kact];
    reCoeff(sx_gg_0,   gg_0,   SEpower, Kact-1);
    reCoeff(sx_gg_1,   gg_1,   SEpower, Kact-1);
    reCoeff(sx_gg_2_1, gg_2_1, SEpower, Kact-1);
    reCoeff(sx_gg_2_0, gg_2_0, SEpower, Kact-1);
    for(int j=0; j<Kact; j++) {
      res -= pow(1-z,j) *( 0+gg_0  [j]
			   + gg_1  [j] * xht
			   + gg_2_0[j] * xht*xht
			   + gg_2_1[j] * xht*xht * lHt );
    }
    res *= pow(z,-SEpower);
    // add exact small x
    res += -9*c2*log(z)/z;
    // damping
    //res * pow(1-z, 10*z/(1-z));
    res *= pow(1-z, Kact+1);
    //
    return res;
  }
  //*/

};

