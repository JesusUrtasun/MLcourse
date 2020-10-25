/*--------------------------

The aim of this library is the computation of the Mellin transforms of

Dk(z)  = [ log^k (1-z) / (1-z) ]_+

DMk(z) = [ log^k (log(1/z)) / log(1/z) ]_+

DBk(z) = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+'

where the prime (') in the last denotes that, in fact,
the log sqrt(z) term is outside the plus-distribution.

---------------------------*/


#include "complex.hh"


dcomplex psi0(dcomplex Z);
dcomplex psi_n(int M, dcomplex Z);
dcomplex psi(dcomplex x);
dcomplex psi(int n, dcomplex x);

double GD(int k, double N);
double DG(int k, double N);
double GB(int k, double N);


//////////////////////////////////////////
// Mellin of Dk(z)  = [ log^k (1-z) / (1-z) ]_+
//
double mell_D(int k, double N);
// derivative wrt N
double mell_D_der(int k, double N);
// complex version
dcomplex mell_D(int k, dcomplex N);
//////////////////////////////////////////
//////////////////////////////////////////
// same as before, but without constant
//
double mell_D_noconst(int k, double N);
// complex version
dcomplex mell_D_noconst(int k, dcomplex N);
//////////////////////////////////////////






//////////////////////////////////////////
// Mellin of DMk(z) = [ log^k (log(1/z)) / log(1/z) ]_+
//  (large N limit of mell_D, i.e. the Minimal prescription logs)
//
double mell_DM(int k, double N);
dcomplex mell_DM(int k, dcomplex N);
//////////////////////////////////////////
//////////////////////////////////////////
// Mellin of DMconstk(z) = [ log^k (log(1/z)) / log(1/z) ]_+ + delta(1-z) Gamma^(k+1) /(k+1)
//  (large N limit of mell_D, i.e. the Minimal prescription logs, including constants)
//
double mell_DMconst(int k, double N);
dcomplex mell_DMconst(int k, dcomplex N);
//////////////////////////////////////////






//////////////////////////////////////////
// Mellin of DBk(z) = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+'
// where prime denotes that the delta term is the same as in D_k(z)
//
double mell_DB(int k, double N);
// complex version
dcomplex mell_DB(int k, dcomplex N);
//////////////////////////////////////////
//////////////////////////////////////////
// same as before, but without constant
//
double mell_DB_noconst(int k, double N);
// complex version
dcomplex mell_DB_noconst(int k, dcomplex N);
//////////////////////////////////////////



//////////////////////////////////////////
// Mellin of DBk(z) - DBtildek = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+' - [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+
// i.e. the constant difference between the two possible defs of the distributions
//
double mell_Bdiff(int k);
//////////////////////////////////////////




//////////////////////////////////////////
// Mellin of Dbark(z) = [ log^k [(1-z)/sqrt(z)] / (1-z) ]_+'  (used in DIS)
// where prime denotes that the delta term is the same as in D_k(z)
//
double mell_Dbar(int k, double N);
// complex version
dcomplex mell_Dbar(int k, dcomplex N);
//////////////////////////////////////////
//////////////////////////////////////////
// same as before, but without constant
//
double mell_Dbar_noconst(int k, double N);
// complex version
dcomplex mell_Dbar_noconst(int k, dcomplex N);
//////////////////////////////////////////

double mell_Dbar_der(int k, double N);




//////////////////////////////////////////
// Mellin of Lpk(z) = (1-z)^p log^k (1-z)
//
double mell_L(unsigned int p, unsigned int k, double N);
// complex version
dcomplex mell_L(unsigned int p, unsigned int k, dcomplex N);
//////////////////////////////////////////

double   mell_L_trunc(unsigned int p, unsigned int k, double   N, unsigned int trunc);
dcomplex mell_L_trunc(unsigned int p, unsigned int k, dcomplex N, unsigned int trunc);


double   mell_Lz(unsigned int p, unsigned int k, double   N);
//dcomplex mell_Lz(unsigned int p, unsigned int k, dcomplex N);
double   mell_Lz_trunc(unsigned int p, unsigned int k, double   N, unsigned int trunc);
//dcomplex mell_Lz_trunc(unsigned int p, unsigned int k, dcomplex N, unsigned int trunc);


double DG(int k, double N);


