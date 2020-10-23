      double complex function harlander(x,t)
      implicit none
      double complex x,t,F0B1,F0B2,HPL2,HPL3,HPL4
      double precision zeta3,pi

      pi = 3.1415926535897932385d0
      zeta3 = 1.2020569031595942854d0

      F0B1 = x*(1.d0+x)**2/(1.d0-x)**4*(72d0*HPL4(1, 0, -1, 0, x) 
     &     + 6d0*log(1.d0-x)*log(x)**3-36d0*pi**2/6d0*HPL2(0, 1, x) 
     &     - 36d0* (pi**2/6d0)* log(1.d0-x)* log(x) 
     &     - 108d0 *zeta3* log(1.d0-x) - 64d0* HPL3(0, 0, 1, -x)  
     &     +32d0* HPL2(0, 1, -x)* log(x) - 8d0* (pi**2/6d0)* log(x))
     &     -36d0* x* (5d0 + 5d0* x + 11d0* x**2 + 11d0* x**3)
     &     /(1.d0-x)**5* HPL4(0, 0, 0, 1, -x) -36d0* x*
     &     (5d0 + 5d0* x + 7d0* x**2 + 7d0* x**3)/(1.d0-x)**5
     &     * HPL4(0, 0, 0, 1, x) + 4d0* x *(1d0 + x)
     &     * (23d0 + 41d0* x**2)/(1.d0-x)**5* (HPL3(0, 0, 1, x) + 
     &     HPL3(0, 0, 1, -x))* log(x)
     &     -16d0* x* (1d0 + x + x**2 + x**3)/(1.d0-x)**5* HPL2(0, 1, -x)
     &     *log(x)**2
     &     -2d0* x* (5d0 + 5d0* x + 23d0* x**2 + 23d0* x**3)/(1.d0-x)**5
     &     *HPL2(0, 1, x)* log(x)**2
     &     +x*(5d0 + 5d0* x - 13d0* x**2 - 13d0* x**3)/24d0/(1.d0-x)**5
     &     * log(x)**4
     &     +x* (1d0 + x - 17d0* x**2 - 17d0* x**3)/(1.d0-x)**5
     &     *(pi**2/6d0)* log(x)**2
     &     +2d0* x* (11d0 + 11d0* x - 43d0* x**2 - 43d0* x**3)
     &     /(1.d0-x)**5* zeta3* log(x)
     &     +36d0* x* (1d0 + x - 3d0* x**2 - 3d0* x**3)/(1.d0-x)**5
     &     * (pi**4/90d0)
     &     -2d0* x* (55d0 + 82d0* x + 55d0* x**2)/(1.d0-x)**4
     &     * HPL3(0, 0, 1, x)
     &     +2d0* x* (51d0 + 74d0* x + 51d0* x**2)/(1.d0-x)**4 
     &     *HPL2(0, 1, x)* log(x)
     &     +x*(47d0 + 66d0* x + 47d0* x**2)/(1.d0-x)**4
     &     * log(1.d0-x)* log(x)**2
     &     +x* (6d0 + 59d0* x + 58d0* x**2 + 33d0* x**3)/3d0/(1.d0-x)**5
     &     * log(x)**3
     &     +2d0* x* (31d0 + 34d0* x + 31d0* x**2)/(1.d0-x)**4
     &     * zeta3
     &     +3d0* x* (3d0 + 22d0* x + 3d0* x**2)/2d0/(1.d0-x)**4
     &     * log(x)**2
     &     -24d0* x* (1d0 + x)/(1.d0-x)**3* log(x) - 94d0*x/(1.d0-x)**2

c      F0B2 = (3d0*(t - ((-2d0 + t)*HPL2(0, 0, x))/2d0 
c     &     + ((-1d0 + t)*log(x))/(2d0*(1d0 + sqrt((-1d0 + t)/t))**2
c     &     *sqrt((-1d0 + t)/t)*t*x)))/t**2

      harlander =  F0B1
c      harlander =  F0B1 !+ 4d0/3d0*F0B2


      return
      end function


