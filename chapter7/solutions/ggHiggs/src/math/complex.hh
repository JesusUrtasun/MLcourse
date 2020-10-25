#pragma once

#include <cmath>
#include <iostream>


class dcomplex{
private:
  double r, i;

public:
  inline dcomplex(double re=0.0, double im= 0.0): r(re), i(im) {}
  inline dcomplex(const dcomplex & a): r(a.r), i(a.i){}

  friend inline double real(const dcomplex & a);
  friend inline double imag(const dcomplex & a);
  friend inline double norm(const dcomplex & a);
  friend inline double arg(const dcomplex & a);
  /*   should be in the range  -PI to PI  */
  friend inline double ABS(const dcomplex & a);
  friend inline double abs(const dcomplex & a);

  inline dcomplex & operator = (const dcomplex & a){
    if (r != a.r || i != a.i) {
      r = a.r; i = a.i;
    }
    return *this;
  }
  inline dcomplex  operator - () const {dcomplex A = *this; A.r = -r; 
    A.i = -i; return A;
  }
  inline dcomplex& operator += (double a){ r += a; return *this;}
  inline dcomplex& operator -= (double a){ r -= a; return *this;}
  inline dcomplex& operator *= (double a){ r *= a; i *= a; return *this;}
  inline dcomplex& operator /= (double a){ r /= a; i /= a; return *this;}
  inline dcomplex& operator += (const dcomplex &a) {
    r +=a.r ; i += a.i; return *this;}
  inline dcomplex& operator -= (const dcomplex &a) {
    r -=a.r; i -= a.i; return *this;}
  inline dcomplex& operator *= (const dcomplex &a) { double rr = r;
    r = rr*a.r - i*a.i; i = i*a.r + rr*a.i; return *this;}
  inline dcomplex& operator /= (const dcomplex &a){
    double m = norm(a); double rr = r;
    r = (rr * a.r+ i*a.i)/m; i = (i*a.r - rr*a.i)/m; return *this;}
  inline dcomplex& operator ++ () {r += 1.0; return *this;}
  inline dcomplex& operator -- () {r -= 1.0; return *this;}

  friend std::istream& operator >> (std::istream& is, dcomplex& a);
  friend std::ostream& operator << (std::ostream& os,  const dcomplex& a);

  friend inline dcomplex operator + (const dcomplex &a, const dcomplex &b);
  friend inline dcomplex operator + (const dcomplex &a, double b);
  friend inline dcomplex operator + (double b, const dcomplex &a);
  friend inline dcomplex operator - (const dcomplex &a, const dcomplex &b);
  friend inline dcomplex operator - (const dcomplex &a, double b);
  friend inline dcomplex operator - (double b, const dcomplex &a);
  friend inline dcomplex operator * (const dcomplex &a, const dcomplex &b);
  friend inline dcomplex operator * (const dcomplex &a, double b);
  friend inline dcomplex operator * (double b, const dcomplex &a);
  friend inline dcomplex operator / (const dcomplex &a, const dcomplex &b);
  friend inline dcomplex operator / (const dcomplex &a, double b);
  friend inline dcomplex operator / (double b, const dcomplex &a);
  friend inline int operator == (const dcomplex &a, const dcomplex &b);
  friend inline int operator != (const dcomplex &a, const dcomplex &b);
  friend inline dcomplex conj(const dcomplex & a);
  friend inline dcomplex exp(const dcomplex & a);
  friend inline dcomplex log(const dcomplex &a);
  friend inline dcomplex cos(const dcomplex & a);
  friend inline dcomplex sin(const dcomplex & a);

  /*  constructor from polar coordinates  */
  friend inline dcomplex polar(double m,double ar);
};

const dcomplex I(0.0,1.0);

inline double real(const dcomplex & a){ return a.r;}
inline double imag(const dcomplex & a){ return a.i;}
inline double norm(const dcomplex & a){ return (a.r*a.r + a.i* a.i);}
inline double arg(const dcomplex & a){ return atan2(a.i,a.r); }

inline dcomplex operator + (const dcomplex &a, const dcomplex &b){
   return dcomplex(a.r+b.r,a.i+b.i);
}
inline dcomplex operator + (const dcomplex &a, double b){
   return dcomplex(a.r + b, a.i);
}
inline dcomplex operator + (double b, const dcomplex &a){
   return dcomplex(a.r + b, a.i);
}

inline dcomplex operator - (const dcomplex &a, const dcomplex &b){
   return dcomplex(a.r-b.r,a.i-b.i);
}
inline dcomplex operator - (const dcomplex &a, double b){
   return dcomplex(a.r - b, a.i);
}
inline dcomplex operator - (double b, const dcomplex &a){
   return dcomplex(b -a.r, -a.i);
}

inline dcomplex operator * (const dcomplex &a, const dcomplex &b){
   return dcomplex(a.r*b.r-a.i*b.i,a.i*b.r + a.r*b.i);
}
inline dcomplex operator * (const dcomplex &a, double b){
   return dcomplex(a.r* b, a.i*b);
}
inline dcomplex operator * (double b, const dcomplex &a){
   return dcomplex(b *a.r, b*a.i);
}

inline dcomplex operator / (const dcomplex &a, const dcomplex &b){
   double m = norm(b);
   return dcomplex((a.r*b.r+a.i*b.i)/m,(a.i*b.r- a.r*b.i)/m);
}
inline dcomplex operator / (const dcomplex &a, double b){
   return dcomplex(a.r/b, a.i/b);
}
inline dcomplex operator / (double b, const dcomplex &a){
   double m = norm(a);
   return dcomplex(b *a.r/m, -b*a.i/m);
}

inline int operator == (const dcomplex &a, const dcomplex &b) {
   return((a.r ==b.r) && (a.i == b.i) ? 1 : 0 );}

inline int operator != (const dcomplex &a, const dcomplex &b) {
   return((a.r == b.r) && (a.i == b.i) ? 0 : 1);}

inline dcomplex polar(double m, double argu){
   return  dcomplex(m*cos(argu),m*sin(argu));}

inline dcomplex conj(const dcomplex & a) {
   return dcomplex(a.r,-a.i);}

inline dcomplex sqrt(const dcomplex & a) {
  double n = norm(a);
  double m = std::sqrt(std::sqrt(n));
  return polar(m,0.5*arg(a));}

inline dcomplex exp(const dcomplex & a){
  double m = std::exp(a.r);
   return  polar(m,a.i);}

inline dcomplex log(const dcomplex &a){
  return dcomplex( 0.5 * std::log(norm(a)), arg(a));}

inline dcomplex pow(double a,const dcomplex & b){
   return exp(log(a)*b);}

inline dcomplex pow(const dcomplex & a,double b){
   return exp(log(a)*b);}

inline dcomplex pow(const dcomplex & a,const dcomplex & b){
   return exp(log(a)*b);}

inline dcomplex cos(const dcomplex & a){
  return dcomplex( std::cos(a.r)*std::cosh(a.i), std::sin(a.r)*std::sinh(a.i)) ;}

inline dcomplex sin(const dcomplex & a){
  return dcomplex( std::sin(a.r)*std::cosh(a.i), std::cos(a.r)*std::sinh(a.i)) ;}

inline double ABS(const dcomplex & a){
  return std::sqrt(a.r*a.r + a.i* a.i);}

inline double abs(const dcomplex & a){
  return std::sqrt(a.r*a.r + a.i* a.i);}



