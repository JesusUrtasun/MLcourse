#include <fstream>
#include <vector>
#include <iostream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "finite-mt_small-x_coeffs.hh"
#include <string>

using namespace std;




SmallxCoeffs::SmallxCoeffs(double mH, double mtop, double mbot, double mcha) {
  hsx = new HiggsSmallX(mH,mtop,mbot,mcha);
}
SmallxCoeffs::~SmallxCoeffs() {
  delete hsx;
}


double SmallxCoeffs::eval(double muFrat, map<double,double> &mc) {
  map<double,double>::iterator it;
  it = mc.find(muFrat);
  if(it == mc.end()) {
    hsx->SetFactorizationScale(muFrat);
    std::vector<std::vector<double> > hk = hsx->hCoeffs(3);
    mc10[muFrat] = hk[1][0];
    mc20[muFrat] = hk[2][0];
    mc11[muFrat] = hk[2][1];
    mc30[muFrat] = hk[3][0];
    mc21[muFrat] = hk[3][1];
    // deprecated
    //mc1[muFrat] = 2*mc10[muFrat];
    //mc2[muFrat] = 2*mc20[muFrat] + mc11[muFrat];
    //mc3[muFrat] = 2*(mc30[muFrat] + mc21[muFrat]);
  }
  return mc[muFrat];
}
double SmallxCoeffs::c10(double muFrat) {
  return eval(muFrat, mc10);
}
double SmallxCoeffs::c20(double muFrat) {
  return eval(muFrat, mc20);
}
double SmallxCoeffs::c11(double muFrat) {
  return eval(muFrat, mc11);
}
double SmallxCoeffs::c30(double muFrat) {
  return eval(muFrat, mc30);
}
double SmallxCoeffs::c21(double muFrat) {
  return eval(muFrat, mc21);
}

/*
// deprecated
double SmallxCoeffs::c1marzani(double muFrat) {
  return eval(muFrat, mc1);
}
double SmallxCoeffs::c2marzani(double muFrat) {
  return eval(muFrat, mc2);
}
double SmallxCoeffs::c3marzani(double muFrat) {
  return eval(muFrat, mc3);
}
*/




