/*
  ggHiggs code based on

    R. D. Ball, M. Bonvini, S. Forte, S. Marzani, G. Ridolfi
    arXiv:1303.3590

  please cite it if you use it
*/

#include <iomanip>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sys/stat.h> 
#include <sys/time.h> 
#include <getopt.h> 

#include "./src/ggHiggs.hh"

using namespace std;
using namespace ggHiggs;

/*
  If you like, you can customize your alpha_s routine.
  However, consistency with PDF sets suggests the use of the LHAPDF alphas routine.
*/
double ggHiggs::alphas(double mu) {
  //return alphasPDF(mu);
  return alphas_ihixs(mu, 3, alphasPDF(GetZmass()), GetZmass()); // taken from ihixs code
  //return alphas_ihixs(mu, order, alphasPDF(GetZmass()), GetZmass()); // taken from ihixs code
}

void print_usage() {
  cout << endl
       << "Usage:" << endl
       << " -p  --PDF-set <set>    specify PDFs being used (default \"PDF4LHC15_nnlo_100\")" << endl
       << " -0  --PDF-member <i>   specify PDF member (default 0)" << endl
       << " -m  --mass <mass>      set Higgs mass [GeV]" << endl
       << " -E  --Ecm <energy>     set collider energy [GeV]" << endl
       << " -f  --fact <murat>     set muF/mH ratio" << endl
       << " -r  --ren <murat>      set muR/mH ratio" << endl
       << " -t  --topmass <mass>   set top quark mass [GeV] (either mt_pole or mt(mt), dependeing on the -S flag)" << endl
       << " -b  --include-bottom   include the bottom Yukawa coupling" << endl
       << " -c  --include-charm    include the charm  Yukawa coupling" << endl
       << " -e  --noEFT            use finite quark masses (default: infinite-mtop EFT approximation)" << endl
#ifdef withHELL
       << " -x  --small-x          include small-x resummation (from HELL-x) in the ggH partonic cross sections" << endl
#endif
       << " -S  --MSbarMass        use MSbar running top mass (with mt(mt) = 162.7 GeV)" << endl
       << " -o  --CPodd            pseudo-scalar (CP odd) Higgs production" << endl
       << " -v  --verbose          show verbose output" << endl
       << endl;
  cout << "Most realistic setup:" << endl
       << "  ggHiggs -e -b -c" << endl
       << endl;
  exit(0);
  return;
}

double muFrat = 1;
double muRrat = 1;
const string sorder[4] = {"LO", "NLO", "NNLO", "NNNLO"};
void read_arguments(int argc, char* argv[], double &mH, double &sqrts, string &set, int &member, bool &nob, bool &noc, bool &smallx) {
  const char* const short_options = "ht:ep:0:f:r:m:E:bcSvox";
  const struct option long_options[] = { { "help", 0, NULL, 'h' },
					 { "topmass", 1, NULL, 't' },
					 { "noEFT", 0, NULL, 'e' },
					 { "PDF-set", 1, NULL, 'p' },
					 { "PDF-member", 1, NULL, '0' },
					 { "fact", 1, NULL, 'f' },
					 { "ren", 1, NULL, 'r' },
					 { "mass", 1, NULL, 'm' },
					 { "Ecm", 1, NULL, 'E' },
					 { "include-bottom", 0, NULL, 'b' },
					 { "include-charm", 0, NULL, 'c' },
					 { "MSbarMass", 0, NULL, 'S' },
					 { "CPodd", 0, NULL, 'o' },
					 { "small-x", 0, NULL, 'x' },
					 { "verbose", 0, NULL, 'v' },
					 { NULL, 0, NULL, 0 } };
  int next_option;
  ostringstream sset;
  do {
    next_option = getopt_long (argc, argv, short_options, long_options, NULL);
    switch (next_option) {
    case 'h':
      print_usage();
    case 'e':
      finite_mt = true;
      break;
    case 't':
      mtop = strtod(optarg, NULL);
      mtmt = mtop;
      break;
    case 'p':
      sset << optarg;
      set = sset.str();
      break;
    case '0':
      member = strtoul(optarg, NULL, 10);
      break;
    case 'f':
      muFrat = strtod(optarg, NULL);
      break;
    case 'r':
      muRrat = strtod(optarg, NULL);
      break;
    case 'm':
      mH = strtod(optarg, NULL);
      break;
    case 'E':
      sqrts = strtod(optarg, NULL);
      break;
    case 'b':
      nob = false;
      break;
    case 'c':
      noc = false;
      break;
    case 'S':
      MSbarMass = true;
      cout << "\033[0;31m WARNING: MSbar scheme for the top mass is a new feature, implemented only for the top in the large-mt EFT at the moment. Use with care.\033[0m" << endl;
      break;
    case 'o':
      SetCPcharge(CPodd);
      break;
    case 'x':
      smallx = true;
      break;
    case 'v':
      quiet = false;
      break;
    case '?':
      print_usage();
    case -1: break;
    default: abort();
    }
  }
  while (next_option != -1);
  //
  return;
}

int main (int argc, char* argv[]) {

  // relative precision for integrals
  INTEGRAL_PRECISION = 2.e-3;

  // default values (modified by runtime arguments)
  double mH = 125;
  double sqrts = 13000;
  //string set = "PDF4LHC15_nnlo_100";
  string set = "NNPDF31_nlo_as_0118";
  int member = 0;

  // read runtime arguments
  bool nob=true, noc=true; // when true, bottom and charm do NOT run into the loop coupling to the Higgs
  bool smallxRes = false;
  read_arguments(argc, argv, mH, sqrts, set, member, nob, noc, smallxRes);
  _muF_mu_ratio  = muFrat;
  _muR_muF_ratio = muRrat/muFrat;
#ifdef withHELL
  if(smallxRes && !finite_mt) {
    cout << "\033[0;31m  Warning: Small-x resummation assumes exact theory, so it's not available within the EFT. Try running with -e -x\033[0m  " << endl;
    smallxRes = false;
  }
#endif

  // Start time counter
  struct timeval t0, t1;
  gettimeofday(&t0,NULL);

  // init PDFs
  cout << "Loading PDF set " << set << endl;
  pdf_init(set, member);

  // masses (negative mass means the quark will not contribute to the loop)
  mbot = (nob ? -1 : 1) * getQuarkMass(5);
  mcha = (noc ? -1 : 1) * getQuarkMass(4);
  // init
  double tau = mH*mH/sqrts/sqrts;
  sxc = new SmallxCoeffs(mH,mtop,mbot,mcha);
  if(smallxRes) SetSmallxResummation(true);
  init_consts(mH);
  SetSmallxMode(0,0);

  cout << "Higgs mass:       " << mH << " GeV" << endl
       << "Collider energy:  " << sqrts << " GeV" << endl
       << "tau = mH^2/s =    " << tau << endl
       << "top    mass =     " << mtop << "  running in the loop" << endl
       << "bottom mass =     " << fabs(mbot) << "  " << (nob ? "NOT " : "") << "running in the loop" << endl
       << "charm  mass =     " << fabs(mcha) << "  " << (noc ? "NOT " : "") << "running in the loop" << endl
       << "quark mass effects beyond LO: " << (finite_mt ? "YES (exact theory)" : "NO (EFT)") << endl
       << "Renormalization scale muR = " << _muR_muF_ratio*_muF_mu_ratio << " * mH" << endl
       << "Factorization   scale muF = " << _muF_mu_ratio << " * mH" << endl
       << "alphas(muR) = " << alphas(_muR_muF_ratio*_muF_mu_ratio*mH) << endl;

  /******* RUN EXAMPLE *******/

  // Note also that the same PDF set has to be used at every order.
  double LO, NLO, NNLO, NNNLO;
  double LOpref = xs_LO_prefactor(mH);
  double muF = mH *_muF_mu_ratio;
  double muR = muF*_muR_muF_ratio;
  double alpha_s = alphas(muR);
  // results
  double rexact[3];

  bool ggonly = false;

  // LO
  LO = tau*Lum(tau,muF)*LOpref*pow(alpha_s,2);
  if(!quiet) cout << "sigma0 = " << LO << endl;

  // NLO
  order = 1;
  __all_channels__ = true & !ggonly;
  compute(mH, sqrts, LO, 0, 0, &rexact[0]);
  NLO = LO + rexact[0];

  // NNLO
  order = 2;
  __all_channels__ = true & !ggonly;
  compute(mH, sqrts, LO, NLO, 0, &rexact[1]);
  NNLO = NLO + rexact[1];

  // NNNLO
  order = 3;
  __all_channels__ = true & !ggonly;
  compute(mH, sqrts, LO, NLO, NNLO, &rexact[2]);
  NNNLO = NNLO + rexact[2];

  if(smallxRes) {
    cout << "Computing small x resum" << endl;
    NLO   += computeSxRes(mH, sqrts, LO, 1);
    NNLO  += computeSxRes(mH, sqrts, LO, 2);
    NNNLO += computeSxRes(mH, sqrts, LO, 3);
  }

  /*
  // final results
  double appr1[3], appr2[3], appr_mean[3], ProcError1[3], ProcError2[3];
  double herr = 0.05; // momentum violation (fraction)
  for(int i=2; i<3; i++) {
    appr1[i] = (rsoft1[i]+rhighenergy[i]);
    appr2[i] = (rsoft2[i]+rhighenergy[i]);
    appr_mean[i] = appr2[i];
    ProcError1[i] = fabs(appr1[i]-appr2[i]);   // estimated error from soft uncertainty only
    ProcError2[i] = herr*fabs(momconsunc[i]);  // estimated error from high-energy uncertainties [momentum violation only]
    ProcError2[i] = sqrt(pow(herr*fabs(momconsunc[i]),2)+pow(rhighenergy[i]-rhighenergyLL[i],2));  // estimated error from high-energy uncertainties
  }
  */

  // Open file for storing data
  cout << "Generating data file" << endl;
  fstream file_1, file_2, file_3;

  // Full tau range - from sqrts = 1000 GeV to 15000 GeV
  // file_1.open("sigma_exact_long.txt", ios::app);

  // // Write in the file
  // if (file_1.good()) {
  //   file_1 << mH << " " << sqrts << " " << tau << " " << LO << " " << NLO << " " << NNLO << " " << NNNLO << endl;
  // }
  // else
  //   cout << "Error: file is not good" << endl;

  // // Close file
  // file_1.close();

  // Large tau - for sqrts < 5000 GeV
  if (sqrts <= 5000) {
    file_2.open("sigma_large_tau_new.txt", ios::app);

    // Write in the file
    if (file_2.good()) {
      file_2 << mH << " " << sqrts << " " << tau << " " << LO << " " << NLO << " " << NNLO << " " << NNNLO << endl;
    }
    else
      cout << "Error: file is not good" << endl;

    // Close file
    file_2.close();
  }

  // Low tau - for sqrts > 10000 GeV
  if (sqrts >= 11000) {
    file_3.open("sigma_low_tau_new.txt", ios::app);

    // Write in the file
    if (file_3.good()) {
      file_3 << mH << " " << sqrts << " " << tau << " " << LO << " " << NLO << " " << NNLO << " " << NNNLO << endl;
    }
    else
      cout << "Error: file is not good" << endl;

    // Close file
    file_3.close();
  }

  cout << endl
       << "The following cross-sections [in pb] have been computed with the same PDF set (" << set << ") at every order."
       << endl << endl;
  cout << "LO    = " << LO << endl;
  cout << "NLO   = " << NLO << endl;
  cout << "NNLO  = " << NNLO << endl;
  //cout << "NNNLO = " << NNLO + appr_mean[2] << " +/- " << ProcError1[2] << "(soft) +/- " << ProcError2[2] << "(high-energy)" << endl;
  //cout << "      = " << NNLO + appr_mean[2] << " +/- " << sqrt(pow(ProcError1[2],2)+pow(ProcError2[2],2)) << "(total, sum in quadrature)" << endl;
  //cout << "      = " << NNLO + appr_mean[2] << " +/- " << ProcError1[2]+ProcError2[2] << "(total, linear sum)" << endl;
  cout << "NNNLO = " << NNNLO << endl;
  cout << endl;

  cout << "K-factor NLO   = " << NLO/LO << endl;
  cout << "K-factor NNLO  = " << NNLO/LO << endl;
  cout << "K-factor NNNLO = " << NNNLO/LO << endl;
  cout << endl;

  double EW = 1+deltaEW(mH,mtop);
  cout << "ElectroWeak correction factor = " << EW << endl << endl;

  cout << "LO(exact) / LO(EFT) = " << LOpref/xs_LO_prefactor_EFT() << endl;

  // Stop time counter
  gettimeofday(&t1,NULL);
  cout << endl << "Total time: " << t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001 << "s" << endl;

  return 0;

}