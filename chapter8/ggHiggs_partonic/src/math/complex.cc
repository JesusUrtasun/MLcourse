#include "complex.hh"
#include <iomanip>


// stream operators
std::istream & operator >> (std::istream& is, dcomplex & a){
  return is >> a.r  >> a.i;
}
std::ostream & operator << (std::ostream & os, const dcomplex & a){
  int sign = a.i>=0 ? 1 : -1;
  std::string ssign = sign>0 ? " + i " : " - i ";
  double im = a.i==0. ? 0. : a.i/sign;
  return os << "( "  << std::setw(11) << a.r 
	    << ssign << std::setw(11) << im  << " )" ;
}

