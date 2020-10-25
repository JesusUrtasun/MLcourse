#pragma once

#include <iostream>
#include <vector>
#include "../math/complex.hh"

namespace ggHiggs {


  // negative bottom/charm mass excludes it from the loop

  extern double G2l(double mh, double mt, double mb=-1, double mc=-1); // Harlander's implementation

  extern double Rqqbar(double z, dcomplex mH2, std::vector<double> m);
  extern double Rqg   (double z, dcomplex mH2, std::vector<double> m, double v);
  extern double Rgg   (double z, dcomplex mH2, std::vector<double> m, double v);

};

