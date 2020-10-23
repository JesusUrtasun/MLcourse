#include <gsl/gsl_integration.h>
#include "integration.hh"




double integration::gauss(double f(double,void*), double xmin, double xmax, double prec, double *error, void *par){
  double result, err;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  F.function = f;
  F.params=par;
  //gsl_integration_qag  (&F, xmin, xmax, 0, prec, 10000, 2, w, &result, &err);
  //gsl_integration_qag  (&F, xmin, xmax, 0, prec, 10000, 6, w, &result, &err);
  gsl_integration_qags (&F, xmin, xmax, 0, prec, 10000,    w, &result, &err);
  gsl_integration_workspace_free (w);
  if(error) *error = err;
  return result;
}

double integration::gauss_inf(double f(double,void*), double xmin, double prec, double *error, void *par){
  double result, err;
  gsl_function F;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);
  F.function = f;
  F.params=par;
  gsl_integration_qagiu (&F, xmin, 0, prec, 100000, w, &result, &err);
  gsl_integration_workspace_free (w);
  if(error) *error = err;
  return result;
}
