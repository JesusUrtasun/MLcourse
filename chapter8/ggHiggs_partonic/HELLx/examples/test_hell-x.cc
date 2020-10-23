/* -----------------------------------------

   Test of HELLx

   ----------------------------------------- */
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sys/time.h> 

#include "hell-x.hh"


using namespace std;




int main (int argc, char* argv[]) {

  struct timeval t0, t1;
  gettimeofday(&t0,NULL);

  int nf = 5;
  HELLx::HELLxnf sxD(nf, HELLx::NLL, "./data/");

  int Npoints = 9, Nas = 3;
  double as = 0.215;
  double x;

  double mQ = 1.4/4.5;

  double x_min = 1e-9;
  double x_max = 0.1;
  for(int j=0; j<Nas; j++) {
    for(int i=0; i<Npoints; i++) {
      x = x_min * exp(i/(Npoints-1.)*log(x_max/x_min));
      HELLx::sqmatrix<double> dP = sxD.DeltaP(as, x, HELLx::NNLO);
      double dKhg = sxD.deltaKhg (as, x, mQ);
      cout << "as = " << as
	   << "     x = "<< setw(6) << x
	   << "    xDeltaPgg = "<< setw(12) << x*dP.gg()
	   << "    xDeltaPqg = "<< setw(12) << x*dP.qg()
	   << "    xDeltaCLg = "<< setw(12) << x*sxD.deltaCLg(as, x, HELLx::NNLO)/nf
	   <<                      setw(12) << x*sxD.deltaMCLg(as, x, mQ)
	   << "    xDeltaC2g = "<< setw(12) << x*sxD.deltaC2g(as, x, HELLx::NNLO)/nf
	   <<                      setw(12) << x*(sxD.deltaMC2g(as, x, mQ) - dKhg)
	   << "    xDeltaKhg = "<< setw(12) << x*dKhg
	   << endl;
    }
    as -= 0.01;
    cout << endl;
  }


  // Finish time
  gettimeofday(&t1,NULL);
  double t=t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001;
  cout << "Total time: " << t << "s" << endl;

  return 0;

}
