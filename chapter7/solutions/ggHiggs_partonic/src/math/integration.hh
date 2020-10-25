#pragma once


namespace integration {

  extern double gauss(double f(double,void*), double xmin, double xmax, double prec, double *error=NULL, void *par=NULL);
  extern double gauss_inf(double f(double,void*), double xmin, double prec, double *error=NULL, void *par=NULL);


};

