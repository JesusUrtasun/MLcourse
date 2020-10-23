#include "bdv.hh"
#include "../parameters.hh"
#include <gsl/gsl_sf.h>


namespace ggHiggs {

// GSL version of complex Li2
dcomplex gslLi2(dcomplex z) {
  gsl_sf_result re, im;
  gsl_sf_complex_dilog_e(abs(z), arg(z), &re, &im);
  return dcomplex(re.val,im.val);
}
// DCOMPLEX version using 't Hooft and Veltman's change of variable
dcomplex CLi2(dcomplex x){
  double x_0 = -0.30;
  double x_1 =  0.25;
  double x_2 =  0.51;
  if (x == 1.) return ZETA2;
  if (real(x) >= x_2) return ZETA2 - CLi2(1.-x) - log(x)*log(1.-x);
  if ((ABS(imag(x)) > 1.) || (real(x)*real(x) + imag(x)*imag(x) > 1.2))
    return - CLi2(1./x) - 0.5 * log(-x) * log(-x) - ZETA2 ;
  if (real(x) <= x_0) {
    dcomplex zz = log(1.-x);
    return -CLi2(-x/(1.-x)) - zz*zz/2. ;
  } else if (real(x) < x_1){
    dcomplex z = - log(1.-x);
    dcomplex temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
                   *(1.-5.*z*z/294.*(1.-7.*z*z/360.
                   *(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
                   *(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)
                    ))))))));
    return temp;
  } else return - CLi2(-x) + CLi2(x*x)/2.;
}
// Li2 selector
dcomplex Li2(dcomplex z) {
  return CLi2(z);
  return gslLi2(z);
}

struct s_complex { double re; double im; };
extern "C" {
  s_complex harlander_(s_complex *x, s_complex *t);
}

extern dcomplex brn_q_ep0(dcomplex x);
double G2l(double mh, double mt, double mb, double mc) {
  std::vector<double> m;
  if(mt>0) m.push_back(mt);
  if(mb>0) m.push_back(mb);
  if(mc>0) m.push_back(mc);
  dcomplex F0H=0, num=0;
  for(int i=0; i<m.size(); i++) {
    dcomplex Wq = 4.*m[i]*m[i]/mh/mh;
    //dcomplex Xq = -Wq/pow(1.+sqrt(1.-Wq),2);
    //F0H += brn_q_ep0(Xq+I*1e-60)*Wq*Wq /64.*(-3./4.);
    F0H += brn_q_ep0(Wq);
    dcomplex t = pow(mh/2./m[i],2) - I*1.e-15;
    dcomplex x = (sqrt(1-1./t)-1)/(sqrt(1-1./t)+1);
    s_complex res, xx={real(x),imag(x)}, tt={real(t),imag(t)};
    res = harlander_(&xx, &tt);
    num += dcomplex(res.re,res.im);
  }
  return real(num/F0H);
}





  // new C++ version  -  Marco
  extern dcomplex sborn_mass(dcomplex mH2, std::vector<double> m);
  dcomplex H00(dcomplex x) {
    dcomplex lx = log(x);
    return lx*lx/2.;
  }
  dcomplex H0(dcomplex x) {
    return log(x);
  }
  dcomplex Li2pol(dcomplex r, dcomplex t) {
    return (Li2(r*exp(I*t)) + Li2(r*exp(-I*t)))/2.;
  }
  dcomplex I3(dcomplex s, dcomplex t, dcomplex u, dcomplex v) {
    if(u==0 || s==0) return 0;
    dcomplex beta = (1+sqrt(1+4*t/u/s))/2;
    if(abs(t/u/s)<1e-5) beta = 1.+t/s/u;
    dcomplex gamm = (1+sqrt(1-4  /v  ))/2;
    dcomplex res = 0;
    if(real(v)<0) {
      res = -Li2(gamm/(gamm+beta-1)) +Li2((gamm-1)/(gamm+beta-1))
	+Li2((beta-gamm)/beta) -Li2((beta-gamm)/(beta-1))
	+(pow(log(beta),2) - pow(log(beta-1),2))/2
	+log(gamm) * log((gamm+beta-1)/beta)
	+log(gamm-1) * log((beta-1)/(gamm+beta-1))
	;
    } else if(real(v)<4.) {
      dcomplex a = sqrt(4/v-1);
      dcomplex r = sqrt((a*a+1)/(a*a+pow(2*beta-1,2)));
      double phi   = acos(real(r*(a*a+2*beta-1)/(1+a*a)));
      double theta = acos(real(r*(a*a-2*beta+1)/(1+a*a)));
      if(abs(beta)-1. < 1e-3) phi = acos(1.-real(2*pow(a*(beta-1.)/(1+a*a),2))); // expansion of the argument of acos
      if(phi<0 || phi>M_PI) exit(55);
      if(theta<0 || theta>M_PI) exit(56);
      res = 2*Li2pol(r,theta) - 2*Li2pol(r,phi) + (phi-theta)*(phi+theta-M_PI);
      //if(res!=res) { std::cout << "\033[0;31m I3 failed for phi=" << phi << "  theta=" << theta << "  a=" << a << "  r=" << r << "\033[0m" << std::endl;}
    } else { // v>4
      res = -Li2(gamm/(gamm+beta-1.)) +Li2((gamm-1.)/(gamm+beta-1.))
	+Li2(gamm/(gamm-beta)) -Li2((gamm-1.)/(gamm-beta))
	+log(gamm/(1.-gamm)) * log((gamm+beta-1.)/(beta-gamm))
	-I*M_PI*log((gamm+beta-1.)/(beta-gamm))
	;
    }
    //if(res!=res) { std::cout << "\033[0;31m I3 failed for s=" << s << "  t=" << t << "  u=" << u << "  v=" << v << "\033[0m" << std::endl;}
    return res *2/(2*beta-1);
  }
  dcomplex W3(dcomplex s, dcomplex t, dcomplex u, dcomplex v) {
    return I3(s,t,u,v) -I3(s,t,u,s) -I3(s,t,u,u);
  }
  dcomplex H3(dcomplex a, dcomplex b, dcomplex c) {
    return -W3(b,a,c,a+b+c);
    //return BDV_H3(a,b,c);
  }
  dcomplex Di(dcomplex s, dcomplex t, dcomplex u, dcomplex y) {
    dcomplex x  = (sqrt(1-4*y)-1)/(sqrt(1-4*y)+1);
    dcomplex xs = (sqrt(1-4/s)-1)/(sqrt(1-4/s)+1);
    dcomplex res = 4 + 4*s/(t+u)*( sqrt(1-4*y)*H0(x) - sqrt(1-4/s)*H0(xs) ) + 8./(t+u)*( H00(x)-H00(xs) );
    return res - 2*( H00(x)-H00(xs) );
  }
  dcomplex Bi(dcomplex s, dcomplex t, dcomplex u, dcomplex y) {
    dcomplex x  = (sqrt(1-4*y)-1)/(sqrt(1-4*y)+1);
    dcomplex xs = (sqrt(1-4/s)-1)/(sqrt(1-4/s)+1);
    dcomplex xt = (sqrt(1-4/t)-1)/(sqrt(1-4/t)+1);
    if(t==0) xt = 1.;
    dcomplex H3sut = H3(s,u,t);
    dcomplex H3tsu = H3(t,s,u);
    dcomplex res = s*(t-s)/(s+t) + 2*t*u*(u+2*s)/pow(s+u,2)*( sqrt(1-4*y)*H0(x) - sqrt(1-4/t)*H0(xt) )
      -(1+t*u/s)*H00(x) + H00(xs) -2*(2*s*s/pow(s+u,2)-1-t*u/s)*(H00(x)-H00(xt))
      +1./2.*(t*u/s+3)*H3sut - H3tsu
      + s/4*(H00(x)-H00(xs)) -(s/2-s*s/(s+u))*(H00(x)-H00(xt)) -s/8*H3sut +s/4*H3tsu;
    //if(res!=res) { std::cout << "\033[0;31m Bi failed for s=" << s << "  t=" << t << "  u=" << u << "  y=" << y << "\033[0m" << std::endl;}
    return res;
  }
  dcomplex Ci(dcomplex s, dcomplex t, dcomplex u, dcomplex y) {
    dcomplex x  = (sqrt(1-4*y)-1)/(sqrt(1-4*y)+1);
    dcomplex xs = (sqrt(1-4/s)-1)/(sqrt(1-4/s)+1);
    dcomplex res = -2*s -H3(u,s,t) + (1./2./y-2)*(H00(x)-H00(xs)) +1./4./y*H3(s,u,t);
    //if(res!=res) { std::cout << "\033[0;31m Ci failed for s=" << s << "  t=" << t << "  u=" << u << "  y=" << y << "\033[0m" << std::endl;}
    return res;
  }
  dcomplex Aqq(dcomplex s, dcomplex t, dcomplex u, dcomplex mH2, std::vector<double> m) {
    dcomplex res = 0;
    for(int i=0; i<m.size(); i++) {
      double m2 = m[i]*m[i];
      dcomplex y = m2/mH2;
      dcomplex si = s/m2;
      dcomplex ti = t/m2;
      dcomplex ui = u/m2;
      res += y * Di(si,ti,ui,y);
    }
    return res/2.;
  }
  dcomplex A2(dcomplex s, dcomplex t, dcomplex u, dcomplex mH2, std::vector<double> m) {
    dcomplex res = 0;
    for(int i=0; i<m.size(); i++) {
      double m2 = m[i]*m[i];
      dcomplex y = m2/mH2;
      dcomplex si = s/m2;
      dcomplex ti = t/m2;
      dcomplex ui = u/m2;
      res += y*y * ( Bi(si,ti,ui,y) + Bi(si,ui,ti,y) );
    }
    //std::cout << res/2. << "  <->  " << BDV_A2(s,t,u) << std::endl;
    return res/2.;
  }
  dcomplex A4(dcomplex s, dcomplex t, dcomplex u, dcomplex mH2, std::vector<double> m) {
    dcomplex res = 0;
    for(int i=0; i<m.size(); i++) {
      double m2 = m[i]*m[i];
      dcomplex y = m2/mH2;
      dcomplex si = s/m2;
      dcomplex ti = t/m2;
      dcomplex ui = u/m2;
      res += y*y * ( Ci(si,ti,ui,y) + Ci(ti,ui,si,y) + Ci(ui,si,ti,y) );
    }
    //std::cout << res/2. << "  <->  " << BDV_A4(s,t,u) << std::endl;
    return res/2.;
  }
  /*
  // this is Spira's implementation, and gives identical results
  dcomplex gS(double t) {
    dcomplex res = 0;
    if(t<=1) res = sqrt(1./t-1.) * asin(sqrt(t));
    else     res = sqrt(1.-1./t) * (log((1.+sqrt(1.-1./t))/(1.-sqrt(1.-1./t)))-I*M_PI) /2;
    return res;
  }
  dcomplex fS(double t) {
    dcomplex res = 0;
    if(t<=1) res = pow(asin(sqrt(t)),2.);
    else     res = -pow(log((1.+sqrt(1.-1./t))/(1.-sqrt(1.-1./t)))-I*M_PI,2.)/4.;
    return res;
  }
  dcomplex Aqqg(dcomplex s, dcomplex mH2, std::vector<double> m) {
    dcomplex res = 0;
    for(int i=0; i<m.size(); i++) {
      double m2 = m[i]*m[i];
      double ro = real(mH2)/m2;
      double si = real(s)/m2;
      res += 2./ro * (1.-2.*si/(si-ro)*(gS(si/4.)-gS(ro/4.)) -(1.+4./(si-ro))*(fS(si/4.)-fS(ro/4.)));
    }
    return res;
  }
  */
  // final results (Rqg and Rgg to be integrated over v in 0<v<1)
  double Rqqbar(double z, dcomplex mH2, std::vector<double> m) {
    dcomplex s = mH2/z;
    double v = 0; // the result depends only on t+u=-s(1-z)=mH^2-s which does not depend on v
    dcomplex t = -s*(1-z)*(1-v);
    dcomplex u = -s*(1-z)*v;
    dcomplex born = sborn_mass(mH2,m)*2./3.;
    double res = 128./27. * z*(1-z) * norm(Aqq(s,t,u,mH2,m)) / norm(born);
    //double res = 128./27. * z*(1-z) * norm(Aqqg(s,mH2,m)) / norm(born);
    return res;
  }
  double Rqg(double z, dcomplex mH2, std::vector<double> m, double v) {
    dcomplex s = mH2/z;
    dcomplex t = -s*(1-z)*(1-v);
    dcomplex u = -s*(1-z)*v;
    dcomplex born = sborn_mass(mH2,m)*2./3.;
    double res = ( (1+(1-z)*(1-z)*v*v)/pow(1-(1-z)*v,2) * 2*z * norm(Aqq(t,s,u,mH2,m)) / norm(born) - (1+(1-z)*(1-z))/z/2. ) / (1-v);
    if(res!=res) { res=0; std::cout << "\033[0;31m Rqg failed for z=" << z << "  v=" << v << "\033[0m" << std::endl;}
    return res + z/2.;
  }
  double Rgg(double z, dcomplex mH2, std::vector<double> m, double v) {
    dcomplex s = mH2/z;
    dcomplex t = -s*(1-z)*(1-v);
    dcomplex u = -s*(1-z)*v;
    dcomplex born = sborn_mass(mH2,m)*2./3.;
    double normAgg = norm(A2(s,t,u,mH2,m)) + norm(A2(u,s,t,mH2,m)) + norm(A2(t,u,s,mH2,m)) + norm(A4(s,t,u,mH2,m));
    double res = ( 8*z*z*z*z * normAgg / norm(born) - pow(1-z+z*z,2) ) / v/(1-v);
    if(res!=res) { res=0; std::cout << "\033[0;31m Rgg failed for z=" << z << "  v=" << v << "\033[0m" << std::endl; }
    return res/z/(1-z);
  }


};
