#include <iostream>
#include "chebyshev.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;



int **Tnk;
unsigned int current_Nmax=0;

// coefficients T_nk of Chebyshev expansion
int T(int n, int k) {
  int value=1;
  if(k>n || n<0 || k<0) {
    value = 0;
    //cerr << "Wrong values of (n,k) in T_{n,k}; possible problems are: n<0, k<0, k>n." << endl;
    //exit(0);
  } else if(((n-k)%2)==1) {
    value = 0;
  } else if(n>0) {
    if(k==0) {
      value= -Tnk[n-2][k];
    } else if(k==n && n!=1) {
      value = 2*Tnk[n-1][k-1];
    } else if(k<n) {
      value = 2*Tnk[n-1][k-1]-Tnk[n-2][k];
    }
  }
  return value;
}

// this function computes coefficients T[n][k] up to n=N
int initT(unsigned int N) {
  if(N<=current_Nmax) return 1;
  else {
    current_Nmax = N;
    Tnk = new int*[N+1];
    for(unsigned int n=0; n<=N; n++) {
      Tnk[n]= new int[N+1];
      for(unsigned int k=0; k<=N; k++) {
	Tnk[n][k]=T(n,k);
      }
    }
    return 1;
  }
  return 0;
}


double ctilde(double *c, int k, int N) {
  double value=0;
  for(int n=k; n<=N; n++){
    value += c[n]*Tnk[n][k];
  }
  if(k==0) value -= 0.5*c[0];
  //if(k==N) value -= 0.5*c[N];  // TEST!!!!
  return value;
}
double cbar(double *c, int p, int N, double M) {
  double value=0;
  for(int k=p; k<=N; k++){
    value += (c[k]*gsl_sf_fact(k))/gsl_sf_fact(k-p);
  }
  return value*pow(-2./M,p);
}



struct sd2 {
  double beta, gamma; int channel;
};


extern double ChebLum(double z, int LumChannel);
double Lum_func(double x, void *p) {
  double z = exp(x);
  sd2 pars = *(sd2*)p;
  return ChebLum(z, pars.channel) * pow(1-z,pars.gamma) * pow(z,pars.beta);
}


// routine adapted from GSL library
int cheb_approx(unsigned int N, const double a, const double b, double *c, void *pars) {
  if(a >= b) {
    cout << "ERROR: in Chebyshev approx, wrong interval [a,b]" << endl;
    abort();
  }
  //
  double bma = 0.5*(b-a);
  double bpa = 0.5*(b+a);
  double fac = 2./(N+1.);
  double *f = new double[N+1];
  //
  for(unsigned int k=0; k<=N; k++) {
    double y = cos(M_PI*(k+0.5)/(N+1.));
    f[k] = Lum_func(y*bma + bpa, pars);
  }
  //
  for(unsigned int j=0; j<=N; j++) {
    double sum = 0.;
    for(unsigned int k=0; k<=N; k++) 
      sum += f[k]*cos(M_PI*j*(k+0.5)/(N+1.));
    c[j] = fac * sum;
  }
  return 1;
}



int ApproxLuminosity(double *coeff, double z_min, unsigned int N, double beta, int gamma, int LumChannel) {
  initT(N);
  sd2 pars = {beta, gamma, LumChannel};
  //
  double M=-log(z_min);
  double *coeff_tmp   = new double[N+1];
  double *coeff_tilde = new double[N+1];
  cheb_approx(N, -M, 0., coeff_tmp, &pars);
  //
  for(unsigned int k=0; k<=N; k++) {
    coeff_tilde[k] = ctilde(coeff_tmp,k,N);
  }
  for(unsigned int k=0; k<=N; k++) {
    coeff[k] = cbar(coeff_tilde,k,N,M);
  }
  return 1;
}
