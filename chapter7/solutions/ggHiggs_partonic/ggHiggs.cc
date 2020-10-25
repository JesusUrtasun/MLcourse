/*
  ggHiggs code based on

    R. D. Ball, M. Bonvini, S. Forte, S. Marzani, G. Ridolfi
    arXiv:1303.3590

  please cite it if you use it
*/

/*
  For reproducing figures in arXiv:1303.3590,
  factorization and renormalization scale need to be set to
  muf = mur = mh, by running <<./ggHiggs -f 1 -r 1>>
  or changing the default values
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
#include <gsl/gsl_sf.h> 
#include <gsl/gsl_integration.h>
#include <cuba.h>

#include "./src/ggHiggs.hh"
#include "./src/math/integration.hh"

using namespace std;
using namespace ggHiggs;

namespace ggHiggs {
  void init_compute(double);
  void initGammaCoeffs(bool gamma_NLL, int nf);
  //
  double Dk(int k, double z);
  extern double _integration_v_;
  double FO_1_delta_const();
  double FO_1_log(double z, double D(int,double));
  double FO_1_nonlog(double z, double mu);
  double FO_1_qg (double z, double mu);
  double FO_1_qqb(double z, double mu);
  //
  double FO_2_delta_const(double mu);
  double FO_2_log(double z, double D(int,double));
  double FO_2_nonlog(double z, double mu);
  double FO_2_qg (double z, double mu);
  double FO_2_qqb(double z, double mu);
  double FO_2_qq (double z, double mu);
  double FO_2_qqp(double z, double mu);
  //
  double FO_3_delta_const(double mu);
  double FO_3_log(double z, double D(int,double));
  double FO_3_nonlog(double z, double mu);
  double FO_3_qg (double z, double mu);
  double FO_3_qqb(double z, double mu);
  double FO_3_qq (double z, double mu);
  double FO_3_qqp(double z, double mu);
  //
  double SxRES_gg(double,double,double,int);
  double SxRES_qg(double,double,double,int);
  double SxRES_qq(double,double,double,int);
  //
  double FO_2_smallx(double z, string chan, int mode, int truncation, int nlogs);
  double FO_2_hard(double, double);
  //double FO_2_gg_smallx(double, double);
  //double FO_2_qg_smallx(double, double);
  //double FO_2_qq_smallx(double, double);
  //
  double FO_3_smallx(double z, string chan, int mode, int truncation, int nlogs);
  double FO_3_hard(double, double);
  //double FO_3_gg_smallx(double, double);
  //double FO_3_qg_smallx(double, double);
  //double FO_3_qq_smallx(double, double);
  //
  //extern double lumi_qqb(double x1, double x2, double mu);
  //extern double lumi_qq (double x1, double x2, double mu);
  //extern double lumi_qqp(double x1, double x2, double mu);
  //extern double lumi_SS (double x1, double x2, double mu);
};

/*
  If you like, you can customize your alpha_s routine.
  However, consistency with PDF sets suggests the use of the LHAPDF alphas routine.
*/
double ggHiggs::alphas(double mu) {
  //return alphasPDF(mu);
  return alphas_ihixs(mu, 3, alphasPDF(GetZmass()), GetZmass()); // taken from ihixs code
  //return alphas_ihixs(mu, order, alphasPDF(GetZmass()), GetZmass()); // taken from ihixs code
}

double mH = 125;
int SXmode = 0;
double api, as;

struct intpars { double N; string name; };
double integrand(double x, void *p) {
  intpars par = *(intpars*)p;
  double N = par.N;
  string name = par.name;
  // change of variables
  double z = x, jac = 1/x;
  //double z = exp(-x/(1-x)), jac = pow(1-x,-2);
  //
  double res = 0;
  if(name == "ggNLO" ) res = api* FO_1_nonlog(z, mH);
  if(name == "ggNNLO") res = api*api* FO_2_nonlog(z, mH);
  if(name == "ggN3LO") res = api*api*api* FO_3_nonlog(z, mH);
  //
  if(name == "ggLLx1") res = SxRES_gg(as, z, _muF_mu_ratio, 1);
  if(name == "ggLLx2") res = SxRES_gg(as, z, _muF_mu_ratio, 2);
  if(name == "ggLLx3") res = SxRES_gg(as, z, _muF_mu_ratio, 3);
  //
  if(name == "ggNNLOunc") res = api*api* ( FO_2_smallx(x, "gg", SXmode, 14, 2) - FO_2_smallx(x, "gg", SXmode, 14, 1) );
  if(name == "ggN3LOunc") res = api*api*api* ( FO_3_smallx(x, "gg", SXmode, 37, 3) - FO_3_smallx(x, "gg", SXmode, 37, 2) );
  //
  return res * pow(z,N) * jac;
}
double Mellin(double N, string name) {
  intpars p = { N, name };
  double res = 0;
  double prec = 2e-4;
  gsl_error_handler_t *old_handler = gsl_set_error_handler (NULL);
  gsl_set_error_handler_off();
  res = integration::gauss(integrand, 0, 1, prec, NULL, &p);
  gsl_set_error_handler (old_handler);
  return res;
}

double integrandPlus(double x, void *p) {
  intpars par = *(intpars*)p;
  double N = par.N;
  string name = par.name;
  // change of variables
  double z = x, jac = 1/x;
  //double z = exp(-x/(1-x)), jac = pow(1-x,-2);
  //
  double res = 0;
  if(name == "ggNLO" ) res = api* FO_1_log(z, Dk);
  if(name == "ggNNLO") res = api*api* FO_2_log(z, Dk);
  if(name == "ggN3LO") res = api*api*api* FO_3_log(z, Dk);
  //
  return res * ( pow(z,N) * jac - 1 );
}
double MellinPlus(double N, string name) {
  intpars p = { N, name };
  double res = 0;
  double prec = 2e-3;
  gsl_error_handler_t *old_handler = gsl_set_error_handler (NULL);
  gsl_set_error_handler_off();
  res = integration::gauss(integrandPlus, 0, 1, prec, NULL, &p);
  gsl_set_error_handler (old_handler);
  //
  if(name == "ggNLO" ) res += api* FO_1_delta_const();
  if(name == "ggNNLO") res += api*api* FO_2_delta_const(mH);
  if(name == "ggN3LO") res += api*api*api* FO_3_delta_const(mH);
  //
  return res;
}

int integrand2d(const int *ndim, const double x[], const int *ncomp, double result[], void *p) {
  intpars par = *(intpars*)p;
  double N = par.N;
  string name = par.name;
  _integration_v_ = x[1];
  // change of variables
  double z = x[0], jac = 1/x[0];
  //double z = exp(x[0]/(1-x[0])), jac = pow(1-x[0],-2);
  //
  double res = 0;
  if(name == "ggNLO" ) res = api* FO_1_nonlog(z, mH);
  //if(name == "ggNNLO") res = api*api* FO_2_nonlog(z, mH);
  //if(name == "ggN3LO") res = api*api*api* FO_3_nonlog(z, mH);
  //if(name == "ggLLx1") res = SxRES_gg(as, z, _muF_mu_ratio, 1);
  //if(name == "ggLLx2") res = SxRES_gg(as, z, _muF_mu_ratio, 2);
  //if(name == "ggLLx3") res = SxRES_gg(as, z, _muF_mu_ratio, 3);
  //
  result[0] = res * pow(z,N) * jac;
  return 0;
}
double Mell2d(double N, string name) {
  intpars p = { N, name };
  double cres[1], error[1], prob[1];
  double epsrel=1e-2, epsabs=1.e-10;
  int last = 4;
  int verbose = 0;
  int nregions, neval, fail;
  //system("export CUBACORES=8");
  int dim = 2;
  Cuhre(dim, 1, integrand2d, &p, 1,
	epsrel, epsabs, verbose | last,
	0, 500000, 9,
	NULL, NULL,
	&nregions, &neval, &fail, cres, error, prob);
  //  
  if(fail == -1) {
    cout << "\033[0;31m  ERROR: dimension out of range.\033[0m" << endl;
  } else if(fail>0) {
    cout << "\033[0;31m  ERROR: precision not reached.\033[0m" << endl;
  }
  return cres[0];
}

void print_usage() {
  cout << endl
       << "Usage:" << endl
       << " -p  --pdf-set <set>    specify PDFs being used (default \"PDF4LHC15_nnlo_100\")" << endl
       << " -m  --mass <mass>      set Higgs mass [GeV]" << endl
       << " -E  --Ecm <energy>     set collider energy [GeV]" << endl
       << " -f  --fact <murat>     set muF/mH ratio" << endl
       << " -r  --ren <murat>      set muR/mH ratio" << endl
       << " -t  --topmass <mass>   set top quark mass [GeV] (either mt_pole or mt(mt), dependeing on the -S flag)" << endl
       << " -b  --include-bottom   include the bottom Yukawa coupling" << endl
       << " -c  --include-charm    include the charm  Yukawa coupling" << endl
       << " -e  --noEFT            use finite quark masses (default: infinite-mtop EFT approximation)" << endl
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

bool useLLp = false;
double muFrat = 0.5;
double muRrat = 0.5;
const string sorder[4] = {"LO", "NLO", "NNLO", "NNNLO"};
void read_arguments(int argc, char* argv[], double &mH, double &sqrts, string &set, bool &nob, bool &noc) {
  const char* const short_options = "ht:ep:f:r:m:E:bcSvoDL";
  const struct option long_options[] = { { "help", 0, NULL, 'h' },
					 { "topmass", 1, NULL, 't' },
					 { "noEFT", 0, NULL, 'e' },
					 { "pdf-set", 1, NULL, 'p' },
					 { "fact", 1, NULL, 'f' },
					 { "ren", 1, NULL, 'r' },
					 { "mass", 1, NULL, 'm' },
					 { "Ecm", 1, NULL, 'E' },
					 { "include-bottom", 0, NULL, 'b' },
					 { "include-charm", 0, NULL, 'c' },
					 { "MSbarMass", 0, NULL, 'S' },
					 { "CPodd", 0, NULL, 'o' },
					 { "verbose", 0, NULL, 'v' },
					 { "damping", 0, NULL, 'D' },
					 { "LLp", 0, NULL, 'L' },
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
    case 'v':
      quiet = false;
      break;
    case 'D':
      SXmode = 1;
      break;
    case 'L':
      cout << "Warning: at the moment this flag does nothing!!" << endl;
      useLLp = true;
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
  double sqrts = 13000;
  string set = "PDF4LHC15_nnlo_100";

  // read runtime arguments
  bool nob=true, noc=true; // when true, bottom and charm do NOT run into the loop coupling to the Higgs
  read_arguments(argc, argv, mH, sqrts, set, nob, noc);
  _muF_mu_ratio  = muFrat;
  _muR_muF_ratio = muRrat/muFrat;
  mtop = 173;
  SetSmallxMode(SXmode);

  // Start time counter
  struct timeval t0, t1;
  gettimeofday(&t0,NULL);

  // init PDFs
  cout << "Loading PDF set " << set << endl;
  int subset = 0;
  pdf_init(set, subset);

  //double x1=0.004, x2=0.05;
  //cout << setw(14) << lumi_qqb(x1, x2, mH) + lumi_qq(x1, x2, mH) + lumi_qqp(x1, x2, mH)
  //     << setw(14) << lumi_SS(x1, x2, mH)
  //     << endl;

  // masses (negative mass means the quark will not contribute to the loop)
  mbot = (nob ? -1 : 1) * getQuarkMass(5);
  mcha = (noc ? -1 : 1) * getQuarkMass(4);
  // init
  double tau = mH*mH/sqrts/sqrts;
  sxc = new SmallxCoeffs(mH,mtop,mbot,mcha);

  // make HELL work

  //SetSmallxResummation(true,"../../../../../repositories/hell/hell-xLLp/data/");
  // SetSmallxResummation(true,HELLdata);  // passed as a gcc macro
  SetSmallxResummation(true, "/home/jesus/Desktop/N3PDF/Tools/Bonvini_gg_Higgs/ggHiggs.v4.0/ggHiggs_Nspace/HELLx/data/");  // passed as a gcc macro
  
  init_consts(mH);
  finite_mt = true;
  init_compute(mH);
  //finite_mt = false;
  initGammaCoeffs(true, 5);
  //ComputeMomConsCoeff(mH);
  as = alphas(_muR_muF_ratio*_muF_mu_ratio*mH);
  api = as/M_PI;
  //
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
  //

  const int Npoints = 500;
  double Nmin = 0.05, Nmax = 12;
  ostringstream oname;
  oname << "output/partonicNspace_" << SXmode << ".dat";
  ofstream ofile(oname.str());

  // Generate data files
  cout << "Generating data files" << endl;
  fstream file_1;

  // Iterate over N
  for(int i=0; i<Npoints; i++) {
    double N = 1+Nmin*exp(i/(Npoints-1.)*log(Nmax/Nmin));
    //
    double ggNLO  = Mell2d(N,"ggNLO" ) + MellinPlus(N,"ggNLO" );
    double ggNNLO = Mellin(N,"ggNNLO") + MellinPlus(N,"ggNNLO");
    double ggN3LO = Mellin(N,"ggN3LO") + MellinPlus(N,"ggN3LO");

    // Correct normalization
    ggNLO = ggNLO / as;
    ggNNLO = ggNNLO / (as * as);
    ggN3LO = ggN3LO / (as * as * as);
    
    //
    double ggNNLOunc = Mellin(N,"ggNNLOunc");
    double ggN3LOunc = Mellin(N,"ggN3LOunc");
    //
    SetSmallxResummationMode(false);
    double ggLLx1 = Mellin(N,"ggLLx1");
    double ggLLx2 = Mellin(N,"ggLLx2");
    double ggLLx3 = Mellin(N,"ggLLx3");
    //
    SetSmallxResummationMode(true);
    double ggLLx1LLp = Mellin(N,"ggLLx1");
    double ggLLx2LLp = Mellin(N,"ggLLx2");
    double ggLLx3LLp = Mellin(N,"ggLLx3");
    //
    double ggNLL  = Mellin(N,"ggNLL");
    double ggNNLL = Mellin(N,"ggNNLL");
    double ggN3LL = Mellin(N,"ggN3LL");
    //

    // Save the values of sigma vs N
    file_1.open("../data/ggHiggs_mesh_long_correct3.txt", ios::app);

    // Write in the file
    if (file_1.good()) {
      file_1 << N << " " << ggNLO << " " << ggNNLO << " " << ggN3LO << endl;
    }
    else
      cout << "Error: file is not good" << endl;

    // Close file
    file_1.close();
    
    /*
    order = 1;
    double NLOgg = 0;//api*FO_nonlog(x,mH);
    order = 2;
    double NNLOgg       = api*api*gg_nnlo_mt_reg(x, mH, mtop, _muF_mu_ratio);
    double NNLOgg_match = api*api*smallx_gg_match_nnlo_z(x, 0, mH/mtop);
    //double NNLOgg_sx    = api*api*smallx_gg_match_nnlo_z(x, c2, mH/mtop) - NNLOgg_match;
    double NNLOgg_sxA   = api*api*FO_2_smallx(x, "gg", SXmode, 14, 1);
    double NNLOgg_sx    = api*api*FO_2_smallx(x, "gg", SXmode, 14, 2);
    //double NNLOgg_sx2   = api*api*FO_2_gg_smallx(x, mH);
    double NNLOgg_sx3   = SxRES_gg(as, x, _muF_mu_ratio, 1)-SxRES_gg(as, x, _muF_mu_ratio, 2);
    double NNLOgg_he    = api*api*FO_2_hard(x, mH);
    double NNLOgg_res   = SxRES_gg(as, x, _muF_mu_ratio, 2);
    //
    double NNLOqg       = api*api*qg_nnlo_mt(x, mH, mtop, _muF_mu_ratio);
    //double NNLOqg_sx    = api*api*smallx_gg_match_nnlo_z(x, (c2-c20)*CF/CA, mH/mtop) - NNLOgg_match;
    double NNLOqg_sxA   = api*api*FO_2_smallx(x, "qg", SXmode, 14, 1);
    double NNLOqg_sx    = api*api*FO_2_smallx(x, "qg", SXmode, 14, 2);
    //double NNLOqg_sx2   = api*api*FO_2_qg_smallx(x, mH);
    double NNLOqg_sx3   = SxRES_qg(as, x, _muF_mu_ratio, 1)-SxRES_qg(as, x, _muF_mu_ratio, 2);
    double NNLOqg_res   = SxRES_qg(as, x, _muF_mu_ratio, 2);
    //
    double NNLOqq       = api*api*qqb_nnlo_mt(x, mH, mtop, _muF_mu_ratio);
    //double NNLOqq_sx    = api*api*smallx_gg_match_nnlo_z(x, (c2-2*c20)*pow(CF/CA,2), mH/mtop) - NNLOgg_match;
    double NNLOqq_sxA   = api*api*FO_2_smallx(x, "qq", SXmode, 14, 1);
    double NNLOqq_sx    = api*api*FO_2_smallx(x, "qq", SXmode, 14, 2);
    //double NNLOqq_sx2   = api*api*FO_2_qq_smallx(x, mH);
    double NNLOqq_sx3   = SxRES_qq(as, x, _muF_mu_ratio, 1)-SxRES_qq(as, x, _muF_mu_ratio, 2);
    double NNLOqq_res   = SxRES_qq(as, x, _muF_mu_ratio, 2);
    //
    order = 3;
    double N3LOgg_he    = api*api*api*FO_3_hard(x, mH);
    double N3LOgg_sxA   = api*api*api*FO_3_smallx(x, "gg", SXmode, 37, 1);
    double N3LOgg_sxB   = api*api*api*FO_3_smallx(x, "gg", SXmode, 37, 2);
    double N3LOgg_sxC   = api*api*api*FO_3_smallx(x, "gg", SXmode, 37, 3);
    //double N3LOgg_sx2   = api*api*api*FO_3_gg_smallx(x, mH);
    double N3LOgg_sx3   = SxRES_gg(as, x, _muF_mu_ratio, 2)-SxRES_gg(as, x, _muF_mu_ratio, 3);
    double N3LOgg       = api*api*api* FO_3_nonlog(x, mH) - N3LOgg_sxC;
    double N3LOgg_res   = SxRES_gg(as, x, _muF_mu_ratio, 3);
    //
    double N3LOqg_sxA   = api*api*api*FO_3_smallx(x, "qg", SXmode, 37, 1);
    double N3LOqg_sxB   = api*api*api*FO_3_smallx(x, "qg", SXmode, 37, 2);
    double N3LOqg_sxC   = api*api*api*FO_3_smallx(x, "qg", SXmode, 37, 3);
    //double N3LOqg_sx2   = api*api*api*FO_3_qg_smallx(x, mH);
    double N3LOqg_sx3   = SxRES_qg(as, x, _muF_mu_ratio, 2)-SxRES_qg(as, x, _muF_mu_ratio, 3);
    double N3LOqg       = api*api*api* FO_3_qg(x, mH) - N3LOqg_sxC;
    double N3LOqg_res   = SxRES_qg(as, x, _muF_mu_ratio, 3);
    //
    double N3LOqq_sxA   = api*api*api*FO_3_smallx(x, "qq", SXmode, 37, 1);
    double N3LOqq_sxB   = api*api*api*FO_3_smallx(x, "qq", SXmode, 37, 2);
    double N3LOqq_sxC   = api*api*api*FO_3_smallx(x, "qq", SXmode, 37, 3);
    //double N3LOqq_sx2   = api*api*api*FO_3_qq_smallx(x, mH);
    double N3LOqq_sx3   = SxRES_qq(as, x, _muF_mu_ratio, 2)-SxRES_qq(as, x, _muF_mu_ratio, 3);
    double N3LOqqb      = api*api*api* FO_3_qqb(x, mH) - N3LOqq_sxC;
    double N3LOqq       = api*api*api* FO_3_qq (x, mH) - N3LOqq_sxC;
    double N3LOqqp      = api*api*api* FO_3_qqp(x, mH) - N3LOqq_sxC;
    double N3LOqq_res   = SxRES_qq(as, x, _muF_mu_ratio, 3);
    //
    //cout << setw(13) << (NNLOgg_sx3-NNLOgg_sx2)/NNLOgg_sx2
    //	 << setw(13) << (NNLOqg_sx3-NNLOqg_sx2)/NNLOqg_sx2
    //	 << setw(13) << (NNLOqq_sx3-NNLOqq_sx2)/NNLOqq_sx2
    //	 << endl;
    //cout << setw(13) << (N3LOgg_sx3-N3LOgg_sx2)/N3LOgg_sx2
    //	 << setw(13) << (N3LOqg_sx3-N3LOqg_sx2)/N3LOqg_sx2
    //	 << setw(13) << (N3LOqq_sx3-N3LOqq_sx2)/N3LOqq_sx2
    //	 << endl;
    */
    //
    ofile << setw(12) << N-1
	  << setw(15) << ggNLO
	  << setw(15) << ggNNLO
	  << setw(15) << ggN3LO
	  << setw(15) << ggLLx1  // 5
	  << setw(15) << ggLLx2
	  << setw(15) << ggLLx3
	  << setw(15) << ggNLL   // 8
	  << setw(15) << ggNNLL
	  << setw(15) << ggN3LL
	  << setw(15) << 0       // 11
	  << setw(15) << ggNNLOunc
	  << setw(15) << ggN3LOunc
	  << setw(15) << ggLLx1LLp  // 14
	  << setw(15) << ggLLx2LLp
	  << setw(15) << ggLLx3LLp
	  << endl;
  }

  // Stop time counter
  gettimeofday(&t1,NULL);
  cout << endl << "Total time: " << t1.tv_sec-t0.tv_sec+(t1.tv_usec-t0.tv_usec)*0.000001 << "s" << endl;

  return 0;

}
