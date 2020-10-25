#pragma once
#include <string>
#include <map>
#include "hHiggs_new_var.hh"



class SmallxCoeffs {
private:
  std::map<double,double> mc10, mc20, mc11, mc30, mc21;
  //std::map<double,double> mc1, mc2, mc3;
  double eval(double, std::map<double,double>&);
  HiggsSmallX *hsx;

public:
  SmallxCoeffs(double mH, double mtop, double mbot=-1, double mcha=-1);
  ~SmallxCoeffs();

  double c10(double muFrat);
  //
  double c20(double muFrat);
  double c11(double muFrat);
  //
  double c30(double muFrat);
  double c21(double muFrat);

  // deprecated
  //double c1marzani(double muFrat);
  //double c2marzani(double muFrat);
  //double c3marzani(double muFrat);

};

