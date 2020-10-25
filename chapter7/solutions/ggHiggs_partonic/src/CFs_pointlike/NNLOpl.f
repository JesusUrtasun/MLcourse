c     MARCO: this is taken from Bonciani-Degrassi-Vicini,
c     and is the result from Anastasiou-Melnikov
c     ATTENTION! this is the partonic cross-section, not the coefficient function

      real*8 function c2regAM(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external LI2,LI3,S12
      real *8 LI2,LI3,S12, REGGGscale, Pi
      integer nf
      double precision z, reggg0a, reggg0b, Z2, Z3, muf_mu_rat,
     -     mur_muf_rat, fps
      real *8 A3s,A2s,A1s,A0s,beta0,AK,tf,tr,alfac2,almu2

      Pi = acos(0d0)
      Z2 = 1.6449340668482264366d0
      Z3 = 1.2020569031595942854d0


      reggg0a = 9*(38*z**2-20*z**3+18*z-39*z**4+14d0+7*z**5)*LI3(Z)
     -/(1-z**2)-18*(z**2+z+1d0)**2/(1d0+z)*S12(z**2)
     -+9*(4*z**4+8*z**3+21*z**2+14*z+7d0)/(1d0+z)*S12(-z)
     --9d0/2*(5*z**5-51*z**4-57*z**3+53*z**2+59*z-11d0)/(1-z**2)*S12(z)
     --9d0/2*(8*z**4+8*z**3-3*z**2-2*z-1d0)/(1d0+z)*LI3(-z)
     --9d0/2*(16d0+13*z**5-40*z**3-67*z**4+64*z**2+36*z)/(1-z**2)
     -*LI2(z)*dlog(z)
     -+9d0/2*(2*z**4-15*z**2-10*z-5d0)/(1+z)*LI2(-z)*dlog(z)
     --9d0/4*(59d0+177*z**2-116*z**3+59*z**4-118*z)/(1-z)
     -*dlog(z)*dlog(1-z)**2+27*(3*z**2+2*z+1)/(1+z)*LI2(-z)*dlog(1+z)
     -+9*(6d0-11*z**3+18*z**2-12*z+6*z**4)/(1-z)*dlog(z)**2*dlog(1-z)
     -+9d0/2*(3d0-8*z**3+3*z**4-6*z+9*z**2)/(1-z)*LI2(z)*dlog(1-z)
     --3d0/2*(7*z-7*z**3+4d0+18*z**2-17*z**4+9*z**5)
     -/(1d0-z**2)*dlog(z)**3
     -+9d0/2*(8*z**4+16*z**3+33*z**2+22*z+11d0)/(1+z)*Z2*dlog(1+z)
     --36*(z**2+z+1d0)**2/(1+z)*LI2(z)*dlog(1+z)
     --9d0/4*(4*z**4+8*z**3+27*z**2+18*z+9d0)/(1d0+z)
     -*dlog(1+z)*dlog(z)**2+(-21d0+63d0/2*z**2-18*z+33d0/2*z**3)
     -*dlog(1+z)*dlog(z)+27d0/2*(3*z**2+2*z+1d0)/(1+z)
     -*dlog(1+z)**2*dlog(z)
     --3d0/4*(-280*z**3+143*z**4+394*z-289d0+21*z**2)*LI2(z)/(1-z)
     -+(-21d0+63d0/2*z**2-18*z+33d0/2*z**3)*LI2(-z)
     -+(-2559d0/4*z**3+1079d0/2*z**2-2687d0/4*z+2027d0/4)*dlog(1-z)
     --3d0/8*(374*z**4-389*z+154d0+699*z**2-827*z**3)
     -/(1-z)*dlog(z)**2+(330*z**3-348*z**2+381*z-297d0)*dlog(1-z)**2
     -+3d0/4*(-1180*z**3+641d0-1238*z+1227*z**2+605*z**4)/(1-z)
     -*dlog(z)*dlog(1-z)
     --72*(2d0-z+z**2)*z*dlog(1-z)**3
     --1d0/8*(4318*z**4-6955*z**3+6447*z**2-5611*z+2333d0)
     -/(1-z)*dlog(z)
     -+3d0/4*(495*z**4-886*z**3+564*z**2-200*z+16d0)/(1-z)*Z2
     -+9*(6*z+18*z**2+2d0+10*z**5-6*z**3-19*z**4)/(1d0-z**2)*Z2*dlog(z)
     --9d0/2*(-48*z**3+23*z**4-46*z+3d0+69*z**2)/(1d0-z)*Z2*dlog(1-z)
     -+9d0/2*(-36d0-15*z**4-52*z+19*z**2+13*z**3+33*z**5)/(1-z**2)*Z3
     -+7539/16d0*z**3-24107d0/48*z**2+22879d0/48*z-18157d0/48


      reggg0b = (31d0/6*z+1d0/6+65d0/12*z**2)*S12(z)
     -+(-31d0/12*z**2+1d0/6-17d0/6*z)*LI3(z)
     -+(47d0/12*z**2+25d0/6*z-1d0/6)*LI2(z)*dlog(z)
     -+(-1d0/12*z**2+z/6-1d0/6)*Z2*dlog(1-z)-4*z*(1+z)*Z2*dlog(z)
     -+(-1d0/6*z+1d0/6+z**2/12)*LI2(z)*dlog(1-z)
     -+(1d0/12-z/12+z**2/24)*dlog(1-z)*dlog(z)*dlog((1-z)/z)
     -+5d0/9*z*(1+z)*dlog(z)**3+(-17*z**2/6-7*z/3-1d0/3)*Z3
     -+(-34d0/9*z**3+2*z**2/3-8*z/3+16d0/9)*(dlog(1-z)**2-Z2)
     --2d0/9*(21*z**2+7*z+25*z**4+17d0-61*z**3)/(1-z)*dlog(z)*dlog(1-z)
     -+(785d0/54*z**3-83*z**2/36+49*z/18-461d0/54)*dlog(1-z)
     -+(-351*z**3+117*z**2+68d0+132*z**4+52*z)/72/(1-z)*dlog(z)**2
     -+(227*z**3+68d0+4*z**4-302*z+21*z**2)/36/(1-z)*LI2(1-z)
     -+(333*z**2+2384*z**4-598*z-3041*z**3+1282d0)/216/(1-z)*dlog(z)
     --8887d0/648*z**3+1267d0/432*z**2-497*z/216+12923d0/1296



c     scale depepndent part
c     MARCO: seems to be wrong!!!!!!

      alfac2 = 2*dlog(muf_mu_rat)
      almu2  = 2*dlog(mur_muf_rat*muf_mu_rat)
      tf = -alfac2
      tr = -almu2

c      fps = 0

c      if(tf.ne.0d0) then

      fps=-79*tf - (62*nf*tf)/9d0 - 42*Pi**2*tf +(4*nf*Pi**2*tf)/9d0 - 
     -  (57*tf**2)/2d0 +(16*nf*tf**2)/3d0 -(1089*tr)/8d0 +(33*nf*tr)/4d0
     -   + 99*tf*tr - 6*nf*tf*tr +(147*tf)/(2d0*(-1+z))- 
     -  (66*Pi**2*tf)/(-1+z) +(1031*tf)/(4d0*z)-(86*nf*tf)/(27d0*z)- 
     -  (99*tf**2)/(4d0*z)-(37*nf*tf**2)/(18d0*z) +(363*tr)/(8d0*z)- 
     -  (11*nf*tr)/(4d0*z)-(99*tf*tr)/(2d0*z) +(3*nf*tf*tr)/z +
     -  (12*Pi**2*tf)/((-1+z)*z)- 82*tf*z +(65*nf*tf*z)/9d0 - 
     -  12*Pi**2*tf*z +(4*nf*Pi**2*tf*z)/9d0 -(51*tf**2*z)/4d0 - 
     -  (17*nf*tf**2*z)/6d0 +(1089*tr*z)/8d0 -(33*nf*tr*z)/4d0 - 
     -  (99*tf*tr*z)/2d0 +3*nf*tf*tr*z +(483*tf*z)/(2 - 2*z) +
     -  (60*Pi**2*tf*z)/(-1+z)-(653*tf*z**2)/4d0+(122*nf*tf*z**2)/27d0 - 
     -  9*Pi**2*tf*z**2 +(99*tf**2*z**2)/4d0 +(37*nf*tf**2*z**2)/18d0 - 
     -  (363*tr*z**2)/8d0 +(11*nf*tr*z**2)/4d0 +(99*tf*tr*z**2)/2d0 - 
     -  3*nf*tf*tr*z**2 +(525*tf*z**2)/(2d0*(-1+z))- 
     -  (24*Pi**2*tf*z**2)/(-1+z) +(189*tf*z**3)/(2 - 2*z) +
     -  (12*Pi**2*tf*z**3)/(-1+z)-(3*Pi**2*tf)/(2d0*(1 +z)) +
     -  48*tf*dlog(1-z) +(28*nf*tf*dlog(1-z))/3d0 - 72*tf**2*dlog(1-z) +
     -  198*tr*dlog(1-z)- 12*nf*tr*dlog(1-z)- 
     -  (168*tf*dlog(1-z))/(-1+z)-(198*tf*dlog(1-z))/z - 
     -  (20*nf*tf*dlog(1-z))/(9d0*z) +(36*tf**2*dlog(1-z))/z - 
     -  (99*tr*dlog(1-z))/z +(6*nf*tr*dlog(1-z))/z +
     -  (33*tf*dlog(1-z))/((-1+z)*z)- 6*tf*z*dlog(1-z)- 
     -  (16*nf*tf*z*dlog(1-z))/3d0 +36*tf**2*z*dlog(1-z)- 
     -  99*tr*z*dlog(1-z) +6*nf*tr*z*dlog(1-z) +
     -  (378*tf*z*dlog(1-z))/(-1+z) +90*tf*z**2*dlog(1-z) +
     -  (20*nf*tf*z**2*dlog(1-z))/9d0 - 36*tf**2*z**2*dlog(1-z) +
     -  99*tr*z**2*dlog(1-z)- 6*nf*tr*z**2*dlog(1-z)- 
     -  (384*tf*z**2*dlog(1-z))/(-1+z) +
     -  (141*tf*z**3*dlog(1-z))/(-1+z)- 72*tf*dlog(1-z)**2 +
     -  (216*tf*dlog(1-z)**2)/(-1+z) +(36*tf*dlog(1-z)**2)/z +
     -  36*tf*z*dlog(1-z)**2 -(216*tf*z*dlog(1-z)**2)/(-1+z)- 
     -  36*tf*z**2*dlog(1-z)**2 +(144*tf*z**2*dlog(1-z)**2)/(-1+z)- 
     -  (72*tf*z**3*dlog(1-z)**2)/(-1+z) +
     -  (72*tf*dlog(1-z)**2)/(z - z**2) +
     -  (36*tf -(4*nf*tf)/3d0 +(9*tf)/(2 - 2*z) +(18*tf)/(-1+z) +
     -     (18*tf)/z +54*tf*z -(4*nf*tf*z)/3d0 -(90*tf*z)/(-1+z)- 
     -     9*tf*z**2 +(108*tf*z**2)/(-1+z)-(54*tf*z**3)/(-1+z) +
     -     (9*tf)/(2d0*(1 +z)))*dlog(z)**2 +
     -  dlog(z)*((57*tf)/2d0 -(23*nf*tf)/3d0 +(2*nf*tf**2)/3d0 - 99*tr +
     -     6*nf*tr +(99*tr)/(2 - 2*z) +(132*tf)/(-1+z)- 
     -     (2*nf*tf)/(-1+z) +(18*tf**2)/(-1+z) +(3*nf*tr)/(-1+z) +
     -     (66*tf)/z +(10*nf*tf)/(9d0*z)-(18*tf**2)/z +(99*tr)/(2d0*z)- 
     -     (3*nf*tr)/z +(129*tf*z)/2d0 - nf*tf*z - 54*tf**2*z +
     -     (2*nf*tf**2*z)/3d0 +(99*tr*z)/2d0 - 3*nf*tr*z - 
     -     (234*tf*z)/(-1+z)- 144*tf*z**2 -(2*nf*tf*z**2)/9d0 +
     -     18*tf**2*z**2 -(99*tr*z**2)/2d0 +3*nf*tr*z**2 +
     -     (222*tf*z**2)/(-1+z)-(87*tf*z**3)/(-1+z) +
     -     108*tf*dlog(1-z) +(54*tf*dlog(1-z))/(-1+z)- 
     -     (54*tf*dlog(1-z))/z +(36*tf*dlog(1-z))/((-1+z)*z)- 
     -     54*tf*z*dlog(1-z) +(180*tf*z*dlog(1-z))/(-1+z) +
     -     54*tf*z**2*dlog(1-z)-(216*tf*z**2*dlog(1-z))/(-1+z) +
     -     (108*tf*z**3*dlog(1-z))/(-1+z) +36*tf*dlog(1 +z) +
     -     (18*tf*dlog(1 +z))/z +18*tf*z*dlog(1 +z) +
     -     18*tf*z**2*dlog(1 +z)-(18*tf*dlog(1 +z))/(1 +z))- 
     -  (36*tf*(1 - 4*z +3*z**2 - 2*z**3 +z**4)*LI2(1-z))/
     -   ((-1+z)*z)- 72*tf*(-1+z)*z*LI2((-1+z)/z) +
     -  (36*tf +(18*tf)/z +18*tf*z +18*tf*z**2 -(18*tf)/(1 +z))*
     -   LI2(-z) +(216*tf -(8*nf*tf)/3d0 +(180*tf)/(-1+z) +
     -     (36*tf)/z +108*tf*z -(8*nf*tf*z)/3d0 -(144*tf*z)/(-1+z) +
     -     36*tf*z**2)*LI2(z)

c      endif
c.....Scale dependent part of Ai coefficients

      AK=67D0/6D0-Pi**2/2D0-5*nf/9D0
      beta0=(33D0-2D0*nf)/6D0

      A0s=-(3D0*AK+33D0-6D0*Pi**2)*alfac2+
     -(3D0/2D0*beta0+6D0*beta0)*alfac2**2-
     -3D0/2D0*beta0*almu2*6D0*alfac2


      A1s=-2D0*(33D0-2D0*nf)*alfac2+36D0*alfac2**2
     -+18D0*beta0*almu2


      A2s=-108D0*alfac2

      A3s=0d0

c.....Correction coming from z*Di

       REGGGscale=z*fps-A3s*dlog(1-z)**3-A2s*dlog(1-z)**2
     --A1s*dlog(1-z)-A0s


      c2regAM = reggg0a + nf*reggg0b + REGGGscale
c      write(*,*) z,'    ',reggg0a

      end


      real*8 function c2regAM_qqp(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external LI2,LI3,S12
      real *8 LI2,LI3,S12, REGGGscale, Pi
      integer nf
      double precision z, Z2, Z3, muf_mu_rat, mur_muf_rat, fps
      double precision REGqqp,REGqqp0,REGqqpscale,almu2,alfac2,tr,tf
      Pi = acos(0d0)
      Z2 = 1.6449340668482264366d0
      Z3 = 1.2020569031595942854d0
      REGqqp0=32d0/9*(z+2)**2*(LI3(z)-S12(z))
     --8d0/3*(z+2)**2*dlog(z)*LI2(z)
     --4d0/27*(z+2)**2*dlog(z)**3-8d0/9*(4*z-6d0+z**2)*LI2(z)
     --32d0/9*(z+3)*(1-z)*dlog(1-z)**2
     -+16d0/3*(z+3)*(1-z)*dlog(z)*dlog(1-z)
     -+8d0/9*(z**2+4*z-3d0)*dlog(z)**2+8*Z2/9*(z+2)**2*dlog(z)
     -+4d0/3*(5*z+17d0)*(1-z)*dlog(1-z)
     -+2d0/9*(29*z**2+44*z-59d0)*dlog(z)+(16d0/3-32*z/9-8*z**2/3)*Z2
     --2d0/9*(11*z+105d0)*(1-z)
      alfac2 = 2*dlog(muf_mu_rat)
      almu2  = 2*dlog(mur_muf_rat*muf_mu_rat)
      tf=-alfac2
      tr=-almu2
      REGqqpscale=-8*tf-(32*Pi**2*tf)/27d0+(16*tf**2)/9d0+
     -  (34*tf)/(3d0*z) - 
     -  (32*Pi**2*tf)/(27d0*z)-(8*tf**2)/(3d0*z)-(10*tf*z)/3d0 - 
     -  (8*Pi**2*tf*z)/27d0+(8*tf**2*z)/9d0+(64*tf*dlog(1-z))/9d0 - 
     -  (32*tf*dlog(1 - z))/(3d0*z)+(32*tf*z*dlog(1-z))/9d0 - 
     -  (32*tf*dlog(z))/9d0-(16*tf**2*dlog(z))/9d0+
     -  (16*tf*dlog(z))/(3d0*z) - 
     -  (16*tf**2*dlog(z))/(9d0*z) - (20*tf*z*dlog(z))/9d0 - 
     -  (4*tf**2*z*dlog(z))/9d0 + (16*tf*dlog(z)**2)/9d0 + 
     -  (16*tf*dlog(z)**2)/(9d0*z) + (4*tf*z*dlog(z)**2)/9d0 + 
     -  (64*tf*LI2(z))/9d0 + (64*tf*LI2(z))/(9d0*z) + 
     -  (16*tf*z*LI2(z))/9d0
      REGqqp=REGqqp0+z*REGqqpscale
      RETURN
      END







c     MARCO: this is taken from ggH@nnlo
c     and is the result from Harlander-Kilgore
c     ATTENTION! this is the partonic cross-section, not the coefficient function

      real*8 function c2regHK(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external deltagga, deltaggf
      double precision deltagga, deltaggf, z, muf_mu_rat, mur_muf_rat,
     -     lfh, lfr
      integer nf

      lfh =  2d0*dlog(muf_mu_rat)
      lfr = -2d0*dlog(mur_muf_rat)

      c2regHK = deltagga(z,lfh,lfr) + nf*deltaggf(z,lfh,lfr)

      end



      real*8 function c2regHK_qg(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external deltaqga, deltaqgf
      double precision deltaqga, deltaqgf, z, muf_mu_rat, mur_muf_rat,
     -     lfh, lfr
      integer nf
      lfh =  2d0*dlog(muf_mu_rat)
      lfr = -2d0*dlog(mur_muf_rat)
      c2regHK_qg = deltaqga(z,lfh,lfr) + nf*deltaqgf(z,lfh,lfr)
      end

      real*8 function c2regHK_qqb(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external deltaqqba, deltaqqbf
      double precision deltaqqba, deltaqqbf, z, muf_mu_rat, mur_muf_rat,
     -     lfh, lfr
      integer nf
      lfh =  2d0*dlog(muf_mu_rat)
      lfr = -2d0*dlog(mur_muf_rat)
      c2regHK_qqb = deltaqqba(z,lfh,lfr) + nf*deltaqqbf(z,lfh,lfr)
      end

      real*8 function c2regHK_qq(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external deltaqqa, deltaqqf
      double precision deltaqqa, deltaqqf, z, muf_mu_rat, mur_muf_rat,
     -     lfh, lfr
      integer nf
      lfh =  2d0*dlog(muf_mu_rat)
      lfr = -2d0*dlog(mur_muf_rat)
      c2regHK_qq = deltaqqa(z,lfh,lfr) + nf*deltaqqf(z,lfh,lfr)
      end

      real*8 function c2regHK_qqp(z,nf,muf_mu_rat,mur_muf_rat)
      implicit none
      external deltaqqpa, deltaqqpf
      double precision deltaqqpa, deltaqqpf, z, muf_mu_rat, mur_muf_rat,
     -     lfh, lfr
      integer nf
      lfh =  2d0*dlog(muf_mu_rat)
      lfr = -2d0*dlog(mur_muf_rat)
      c2regHK_qqp = deltaqqpa(z,lfh,lfr) + nf*deltaqqpf(z,lfh,lfr)
      end


      real*8 function deltagga(xx,lfh,lfr)
      implicit none
c      implicit real*8 (a-h,o-z)
      external dilog,trilog
      complex*16 dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f



      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0


      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      deltagga = 67.33333333333333d0 + (1983*dlnx)/4.d0 + 57*dlnx**3 -
     &     (1109*dlxm1)/4.d0 -198*dlnx*dlxm1 + 36*dlnx**2*dlxm1 + 66
     &     *dlxm1**2 - 9*dlnx*dlxm1**2 -144*dlxm1**3 + (21*dlnx**3)/(2
     &     .d0*(2 - xm1)) -(9*dlnx**2*dlxm1)/(2 - xm1) - (139*dlnx)/(2
     &     .d0*xm1) -(33*dlnx**2)/(8.d0*xm1) + (3*dlnx**3)/xm1 + (33
     &     *dlnx*dlxm1)/xm1 +(72*dlnx**2*dlxm1)/xm1 - (279*dlnx*dlxm1**2
     &     )/(2.d0*xm1) -(3435*xm1)/4.d0 - (5823*dlnx*xm1)/4.d0 - (657
     &     *dlnx**2*xm1)/4.d0 -(327*dlnx**3*xm1)/4.d0 + 1530*dlxm1*xm1 +
     &     1017*dlnx*dlxm1*xm1 -45*dlnx**2*dlxm1*xm1 - 675*dlxm1**2*xm1
     &     +(27*dlnx*dlxm1**2*xm1)/2.d0 + 216*dlxm1**3*xm1 + (11107*xm1
     &     **2)/12.d0 +(10365*dlnx*xm1**2)/8.d0 + (2007*dlnx**2*xm1**2)
     &     /8.d0 +24*dlnx**3*xm1**2 - (5567*dlxm1*xm1**2)/4.d0 -1149
     &     *dlnx*dlxm1*xm1**2 - 72*dlnx**2*dlxm1*xm1**2 

      deltagga = deltagga + 642*dlxm1**2*xm1**2 + 135*dlnx*dlxm1**2*xm1
     &     **2 - 144*dlxm1**3*xm1**2 - (7583*xm1**3)/16.d0 - (2171*dlnx
     &     *xm1**3)/4.d0 - (561*dlnx**2*xm1**3)/4.d0 + 3*dlnx**3*xm1**3
     &     + (2583*dlxm1*xm1**3)/4.d0 + 561*dlnx*dlxm1*xm1**3 + 81*dlnx
     &     **2*dlxm1*xm1**3 - 330*dlxm1**2*xm1**3 - (279*dlnx*dlxm1**2
     &     *xm1**3)/2.d0 + 72*dlxm1**3*xm1**3 - 66*z2 + 126*dlnx*z2 +
     &     180*dlxm1*z2 + (81*dlnx*z2)/xm1 + 675*xm1*z2 - 189*dlnx*xm1
     &     *z2 - 270*dlxm1*xm1*z2 - 642*xm1**2*z2 - 18*dlnx*xm1**2*z2 +
     &     180*dlxm1*xm1**2*z2 + (1089*xm1**3*z2)/4.d0 + 81*dlnx*xm1**3
     &     *z2 - 90*dlxm1*xm1**3*z2 - 351*z3 + (1053*xm1*z3)/2.d0 - 351
     &     *xm1**2*z3 + (351*xm1**3*z3)/2.d0 + (-117 + 333*dlnx - 432
     &     *dlxm1 - (27*dlnx)/(2 - xm1) + (36*dlxm1)/(2 - xm1) - 33/(4
     &     .d0*xm1) + (99*dlnx)/(2.d0*xm1) - (261*xm1)/4.d0 - (999*dlnx
     &     *xm1)/2.d0 + 612*dlxm1*xm1 - 138*xm1**2 + 144*dlnx*xm1**2 -
     &     144*dlxm1*xm1**2 + (363*xm1**3)/4.d0 + (99*dlnx*xm1**3)/2.d0
     &     - 36*dlxm1*xm1**3) *dli2a + (-4.5d0 + (45*dlnx)/2.d0 - 72
     &     *dlxm1 + (18*dlnx) /(2 - xm1) - (18*dlxm1)/(2 - xm1) + (189
     &     *xm1)/4.d0 - (99*dlnx *xm1)/4.d0 + 126*dlxm1*xm1 - (81*xm1**2
     &     )/2.d0 - 9*dlnx*xm1**2 - 72*dlxm1*xm1**2 + (33*xm1**3)/4.d0 +
     &     (9*dlnx*xm1**3)/2.d0 + 18*dlxm1*xm1**3)*dli2b 

      deltagga = deltagga + (189 + 27/(2 - xm1) - 171 /(2.d0*xm1) - 234
     &     *xm1 + 36*xm1**2 - (81*xm1**3)/2.d0)*dli3a + (99 + 45 /(2 -
     &     xm1) - (333*xm1)/2.d0 + 72*xm1**2 - 18*xm1 **3)* dli3b + (-99
     &     - 45/(2 - xm1) + (333*xm1) /2.d0 - 72*xm1**2 + 18*xm1 **3)
     &     * dli3c + (-441 - 27/(2 - xm1) - 81/xm1 + (1449*xm1)/2.d0 -
     &     270*xm1**2 - 27*xm1 **3)*dli3d + (29.25d0 - 27/(4.d0*(2 - xm1
     &     )) - (531 *xm1)/8.d0 + 63*xm1**2 - 18*xm1**3)* dli3e + (6
     &     .75d0 - 9/(4.d0*(2 - xm1)) - (189*xm1)/8.d0 + 27*xm1**2 - 9
     &     *xm1**3) *dli3f


c..   lfh-terms:
      deltagga = deltagga + lfh**2*(16.5d0 - 36*dlnx - 72*dlxm1 - (18
     &     *dlnx)/xm1 - (675*xm1)/4.d0 +54*dlnx*xm1 + 108*dlxm1*xm1 +
     &     (321*xm1**2)/2.d0 - 72*dlxm1*xm1**2 -(297*xm1**3)/4.d0 - 18
     &     *dlnx*xm1**3 + 36*dlxm1*xm1**3) +lfh*(139 + 102*dlnx - 45
     &     *dlnx**2 - 66*dlxm1 + 36*dlnx*dlxm1 +216*dlxm1**2 + (9*dlnx
     &     **2)/(2.d0*(2 - xm1)) - (33*dlnx)/(2.d0*xm1) -(45*dlnx**2)/(2
     &     .d0*xm1) + (126*dlnx*dlxm1)/xm1 - (3051*xm1)/4.d0 -513*dlnx
     &     *xm1 + 63*dlnx**2*xm1 + 675*dlxm1*xm1 -54*dlnx*dlxm1*xm1 -
     &     324*dlxm1**2*xm1 + (2773*xm1**2)/4.d0 +(1185*dlnx*xm1**2)/2
     &     .d0 + 9*dlnx**2*xm1**2 - 642*dlxm1*xm1**2 -108*dlnx*dlxm1*xm1
     &     **2 + 216*dlxm1**2*xm1**2 - (2449*xm1**3)/8.d0 -(561*dlnx*xm1
     &     **3)/2.d0 - 27*dlnx**2*xm1**3 + 330*dlxm1*xm1**3 +126*dlnx
     &     *dlxm1*xm1**3 - 108*dlxm1**2*xm1**3 - 90*z2 +135*xm1*z2 - 90
     &     *xm1**2*z2 + 45*xm1**3*z2 +(216 - 18/(2 - xm1) - 306*xm1 + 72
     &     *xm1**2 + 18*xm1**3)*dli2a +(36 + 9/(2 - xm1) - 63*xm1 +
     &     36*xm1**2 - 9*xm1**3)*dli2b)

c..   lfr-terms:
      deltagga = deltagga + lfr*(-99*dlnx + 198*dlxm1 + (99*dlnx)/(2.d0
     &     *xm1) + (297*dlnx*xm1)/2.d0 -297*dlxm1*xm1 - 99*dlnx*xm1**2 +
     &     198*dlxm1*xm1**2 +(363*xm1**3)/8.d0 + (99*dlnx*xm1**3)/2.d0 -
     &     99 *dlxm1*xm1**3 +lfh*(-99 + (297*xm1)/2.d0 - 99*xm1**2 + (99
     &     *xm1 **3)/2.d0))

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltagga = deltagga - ((69*dlnx)/2.d0 - 9*dlnx**2 - 12*dlxm1 +
     &        6*lfh - (3*dlnx)/xm1 +27*xm1 - (39*dlnx*xm1)/2.d0 + 9*dlnx
     &        **2*xm1 + 18*dlxm1*xm1 -9*lfh*xm1 + (57*xm1**2)/4.d0 + 6
     &        *dlnx*xm1**2 - 12*dlxm1*xm1 **2 +6*lfh*xm1**2 - (11*xm1**3
     &        )/4.d0 - 3*dlnx*xm1**3 + 6 *dlxm1*xm1**3 -3*lfh*xm1**3 )

c..   no additional lfr-terms
         
c      endif

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function deltaggf:

      real*8 function deltaggf(xx,lfh,lfr)
      implicit none
c      implicit real*8 (a-h,o-z)
      external dilog,trilog
      complex*16 dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      deltaggf = -3.111111111111111d0 - (265*dlnx)/216.d0 + (287*dlnx**2
     &     )/72.d0 +(17*dlnx**3)/72.d0 + (77*dlxm1)/12.d0 - (68*dlnx
     &     *dlxm1)/9.d0 -(16*dlnx**2*dlxm1)/3.d0 - 4*dlxm1**2 + (16*dlnx
     &     *dlxm1**2)/3.d0 +(5*dlnx)/(3.d0*xm1) + dlnx**2/(4.d0*xm1) -
     &     (2*dlnx*dlxm1)/xm1 +(8441*xm1)/216.d0 + (883*dlnx*xm1)/36.d0
     &     - (8*dlnx**2*xm1)/3.d0 -(dlnx**3*xm1)/3.d0 - (751*dlxm1*xm1)
     &     /18.d0 + (8*dlnx*dlxm1*xm1)/3.d0 +8*dlnx**2*dlxm1*xm1 + (38
     &     *dlxm1**2*xm1)/3.d0 - 8*dlnx*dlxm1**2*xm1 -(17227*xm1**2)/432
     &     .d0 - (2165*dlnx*xm1**2)/72.d0 -(59*dlnx**2*xm1**2)/24.d0 +
     &     (dlnx**3*xm1**2)/8.d0 +(1487*dlxm1*xm1**2)/36.d0 + (26*dlnx
     &     *dlxm1*xm1**2)/3.d0 -(8*dlnx**2*dlxm1*xm1**2)/3.d0 - (32
     &     *dlxm1**2*xm1**2)/3.d0 +(8*dlnx*dlxm1**2*xm1**2)/3.d0 + (8887
     &     *xm1**3)/648.d0 

      deltaggf = deltaggf + (298*dlnx*xm1**3)/27.d0 + (11*dlnx**2*xm1**3
     &     )/6.d0 - (785*dlxm1*xm1**3)/54.d0 - (50*dlnx*dlxm1*xm1**3)/9
     &     .d0 + (34*dlxm1**2*xm1**3)/9.d0 + 4*z2 - (16*dlnx*z2)/3.d0 -
     &     (38*xm1*z2)/3.d0 + 8*dlnx*xm1*z2 + (32*xm1**2*z2)/3.d0 - (8
     &     *dlnx*xm1**2*z2)/3.d0 - (34*xm1**3*z2)/9.d0 + (-12
     &     .13888888888889d0 - (95*dlnx)/12.d0 + (32*dlxm1)/3.d0 + 1/(2
     &     .d0*xm1) + (121*xm1)/6.d0 + 12*dlnx*xm1 - 16*dlxm1*xm1 - (27
     &     *xm1**2)/4.d0 - (47*dlnx*xm1**2)/12.d0 + (16*dlxm1*xm1**2)/3
     &     .d0 + xm1**3/9.d0)*dli2a + (-5.5d0 + 8*xm1 - (17*xm1**2)/6
     &     .d0)*dli3a + (5.25d0 - 8*xm1 + (31*xm1**2)/12.d0)*dli3d


c..   lfh-terms:
      deltaggf = deltaggf + lfh**2*(-1 + (4*dlnx)/3.d0 + (19*xm1)/6.d0 -
     &     2*dlnx*xm1 -(8*xm1**2)/3.d0 + (2*dlnx*xm1**2)/3.d0 + (17*xm1
     &     **3)/18.d0) +lfh*(-3.3333333333333335d0 + (34*dlnx)/9.d0 + (8
     &     *dlnx**2)/3.d0 + 4*dlxm1 -(16*dlnx*dlxm1)/3.d0 + dlnx/xm1 +
     &     (190*xm1)/9.d0 -(4*dlnx*xm1)/3.d0 -4*dlnx**2*xm1 - (38*dlxm1
     &     *xm1)/3.d0 + 8*dlnx*dlxm1*xm1 -(187*xm1**2)/9.d0 - (13*dlnx
     &     *xm1**2)/3.d0 + (4*dlnx**2*xm1**2)/3.d0 +(32*dlxm1*xm1**2)/3
     &     .d0 - (8*dlnx*dlxm1*xm1**2)/3.d0 +(785*xm1**3)/108.d0 + (25
     &     *dlnx*xm1**3)/9.d0 - (34*dlxm1*xm1**3)/9.d0 +(-5
     &     .333333333333333d0 + 8*xm1 - (8*xm1**2)/3.d0)*dli2a)

c..   lfr-terms:
      deltaggf = deltaggf + lfr*(6*dlnx - 12*dlxm1 - (3*dlnx)/xm1 - 9
     &     *dlnx*xm1 + 18*dlxm1*xm1+6*dlnx*xm1**2 - 12*dlxm1*xm1**2 -
     &     (11*xm1**3)/4.d0 - 3*dlnx*xm1**3 +6*dlxm1*xm1**3 + lfh*(6 - 9
     &     *xm1 + 6*xm1**2 - 3*xm1**3))

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaggf = deltaggf - (dlnx + (2*dlnx**2)/3.d0 + (3*xm1)/2.d0 -
     &        dlnx*xm1 - (2*dlnx**2 *xm1)/3.d0 -(5*xm1**2)/3.d0)
c..   no additional lfh- and lfr-terms
c      endif

c-- checked against checks.m [17/06/09,rh]

      end



      real*8 function deltaqga(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      deltaqga = 0.7407407407407407d0 + (6103*dlnx)/108.d0 + (451*dlnx
     &     **2)/54.d0 +(50*dlnx**3)/3.d0 + (353*dlxm1)/18.d0 - (601*dlnx
     &     *dlxm1)/27.d0 +(392*dlnx**2*dlxm1)/9.d0 + (85*dlxm1**2)/36.d0
     &     -(1385*dlnx*dlxm1**2)/18.d0 + (367*dlxm1**3)/54.d0 - (11119
     &     *xm1)/108.d0 -(8521*dlnx*xm1)/54.d0 - (251*dlnx**2*xm1)/9.d0
     &     - 16*dlnx**3*xm1 +(4136*dlxm1*xm1)/27.d0 + (964*dlnx*dlxm1
     &     *xm1)/9.d0 -44*dlnx**2*dlxm1*xm1 - (841*dlxm1**2*xm1)/9.d0
     &     +72*dlnx*dlxm1**2*xm1 + (1153*xm1**2)/72.d0 + (1471*dlnx*xm1
     &     **2)/36.d0 +(74*dlnx**2*xm1**2)/9.d0 + (115*dlnx**3*xm1**2)
     &     /27.d0 -(815*dlxm1*xm1**2)/27.d0 - 38*dlnx*dlxm1*xm1**2 +(52
     &     *dlnx**2*dlxm1*xm1**2)/3.d0 + (325*dlxm1**2*xm1**2)/12.d0
     &     -(553*dlnx*dlxm1**2*xm1**2)/18.d0 + (367*dlxm1**3*xm1**2)/54
     &     .d0 -(1537*xm1**3)/486.d0 - (616*dlnx*xm1**3)/81.d0 -(113
     &     *dlnx**2*xm1**3)/27.d0 + (392*dlxm1*xm1**3)/81.d0 +(400*dlnx
     &     *dlxm1*xm1**3)/27.d0 - 8*dlxm1**2*xm1**3 + (29*z2)/6.d0 

      deltaqga = deltaqga + (629*dlnx*z2)/9.d0 - (50*dlxm1*z2)/9.d0 +
     &     (701*xm1*z2)/9.d0 - 72*dlnx*xm1*z2 - (281*xm1**2*z2)/9.d0 +
     &     (71*dlnx*xm1**2*z2)/3.d0 - (50*dlxm1*xm1**2*z2)/9.d0 + 8*xm1
     &     **3*z2 + (311*z3)/18.d0 + (311*xm1**2*z3)/18.d0 + (-37
     &     .333333333333336d0 + (761*dlnx)/9.d0 - (322*dlxm1)/3.d0 -
     &     (209*xm1)/9.d0 - 96*dlnx*xm1 + 128*dlxm1*xm1 + (59*xm1**2)/18
     &     .d0 + (245*dlnx*xm1**2)/9.d0 - (278*dlxm1*xm1**2)/9.d0 + (26
     &     *xm1**3)/9.d0)*dli2a + (7.87037037037037d0 + 10*dlnx - 10
     &     *dlxm1 - (50*xm1)/9.d0 - 8*dlnx*xm1 + 8*dlxm1*xm1 + (5*xm1**2
     &     )/6.d0 + 2*dlnx*xm1**2 - 2*dlxm1*xm1**2 - (2*xm1**3)/27.d0)
     &     *dli2b + (22.444444444444443d0 - 32*xm1 + (2*xm1**2) /9.d0)
     &     *dli3a + (20 - 16*xm1 + 4*xm1**2)*dli3b + (-20 + 16*xm1 - 4
     &     *xm1**2)*dli3c + (-123 .88888888888889d0 + 128*xm1 - (113*xm1
     &     **2)/3.d0)*dli3d + (-2.5d0 + 2*xm1 - xm1**2/2.d0)*dli3e + (-2
     &     .5d0 + 2*xm1 - xm1**2/2.d0)*dli3f

c..   lfh-terms:
      deltaqga = deltaqga + lfh**2*(-1.5d0 - (160*dlnx)/9.d0 + (31*dlxm1
     &     )/9.d0- (191*xm1)/9.d0 +18*dlnx*xm1 + (50*xm1**2)/9.d0 - (56
     &     *dlnx*xm1**2)/9.d0 +(31*dlxm1*xm1**2)/9.d0 - 2*xm1**3) +lfh*(
     &     -10.444444444444445d0 + (263*dlnx)/27.d0 - (185*dlnx**2)/9.d0
     &     -dlxm1/3.d0 + 76*dlnx*dlxm1 - (31*dlxm1**2)/3.d0 - (2092*xm1)
     &     /27.d0 -(152*dlnx*xm1)/3.d0 + 22*dlnx**2*xm1 + (836*dlxm1*xm1
     &     )/9.d0 -72*dlnx*dlxm1*xm1 + (743*xm1**2)/54.d0 + (52*dlnx*xm1
     &     **2)/3.d0 -(67*dlnx**2*xm1**2)/9.d0 - (218*dlxm1*xm1**2)/9.d0
     &     +(268*dlnx*dlxm1*xm1**2)/9.d0 - (31*dlxm1**2*xm1**2)/3.d0
     &     -(196*xm1**3)/81.d0 - (200*dlnx*xm1**3)/27.d0 + 8*dlxm1*xm1
     &     **3 +(25*z2)/9.d0 + (25*xm1**2*z2)/9.d0 +(54.22222222222222d0
     &     - 64*xm1 + 16*xm1**2)*dli2a+(5 - 4*xm1 + xm1**2)*dli2b)

c..   lfr-terms:
      deltaqga = deltaqga + lfr*(-5.5d0 + (11*dlnx)/2.d0 - 11*dlxm1 + 11
     &     *xm1 + (11*xm1**2)/4.d0 +(11*dlnx*xm1**2)/2.d0 - 11*dlxm1*xm1
     &     **2 + lfh*(5.5d0 + (11*xm1**2)/2.d0))

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaqga = deltaqga - (
     &        0.3333333333333333d0 + 17*dlnx - (28*dlnx**2)/9.d0 + (2
     &        *dlxm1)/3.d0 -lfh/3.d0 + (140*xm1)/9.d0 - (28*dlnx*xm1)/3
     &        .d0 + (28*dlnx**2*xm1)/9.d0 +(17*xm1**2)/6.d0 - (dlnx*xm1
     &        **2)/3.d0 + (2*dlxm1*xm1**2)/3.d0 -(lfh*xm1**2)/3.d0 )
c      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function deltaqgf:

      real*8 function deltaqgf(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

c..   pseudo-scalar:
      deltaqgf = 0.16049382716049382d0 + (10*dlnx)/27.d0 + dlnx**2/9.d0
     &     - (2*dlxm1)/3.d0 -(2*dlnx*dlxm1)/9.d0 + dlxm1**2/18.d0 + (10
     &     *xm1)/27.d0 + (2*dlxm1*xm1)/9.d0 +(179*xm1**2)/162.d0 + (19
     &     *dlnx*xm1**2)/27.d0 + (dlnx**2*xm1**2)/9.d0 -dlxm1*xm1**2 -
     &     (2*dlnx*dlxm1*xm1**2)/9.d0 + (dlxm1**2*xm1**2)/18.d0

c..   lfh-terms:
      deltaqgf = deltaqgf + lfh**2*(0.1111111111111111d0 + xm1**2/9.d0)
     &     +lfh*(0.37037037037037035d0 + (2*dlnx)/9.d0 - (2*dlxm1)/9.d0
     &     + (19*xm1**2)/27.d0 + (2*dlnx*xm1**2)/9.d0 - (2*dlxm1*xm1**2)
     &     /9.d0)

c..   lfr-terms:
      deltaqgf = deltaqgf + lfr*(0.3333333333333333d0 - dlnx/3.d0 + (2
     &     *dlxm1)/3.d0 - (2*xm1)/3.d0 -xm1**2/6.d0 - (dlnx*xm1**2)/3.d0
     &     + (2*dlxm1*xm1**2)/3.d0 +lfh*(-0.3333333333333333d0 - xm1**2
     &     /3.d0))

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function deltaqqa:

      real*8 function deltaqqa(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      deltaqqa = (304*dlnx)/27.d0 + (4*dlnx**2)/27.d0 + (328*dlnx**3)/81
     &     .d0 -(8*dlnx*dlxm1)/9.d0 + 8*dlnx**2*dlxm1 - 16*dlnx*dlxm1**2
     &     -(476*xm1)/27.d0 - (700*dlnx*xm1)/27.d0 - (112*dlnx**2*xm1)
     &     /27.d0 -(8*dlnx**3*xm1)/3.d0 + (88*dlxm1*xm1)/3.d0 + 16*dlnx
     &     *dlxm1*xm1 -(16*dlnx**2*dlxm1*xm1)/3.d0 - (128*dlxm1**2*xm1)
     &     /9.d0 +(32*dlnx*dlxm1**2*xm1)/3.d0 + (44*xm1**2)/9.d0 + (46
     &     *dlnx*xm1**2)/9.d0 +(4*dlnx**2*xm1**2)/3.d0 + (40*dlnx**3*xm1
     &     **2)/81.d0 -(20*dlxm1*xm1**2)/3.d0 - (40*dlnx*dlxm1*xm1**2)/9
     &     .d0 +(8*dlnx**2*dlxm1*xm1**2)/9.d0 + (32*dlxm1**2*xm1**2)/9
     &     .d0 -(16*dlnx*dlxm1**2*xm1**2)/9.d0 + 16*dlnx*z2 + (128*xm1
     &     *z2)/9.d0 -(32*dlnx*xm1*z2)/3.d0 - (32*xm1**2*z2)/9.d0 + (16
     &     *dlnx*xm1**2*z2)/9.d0 +(-0.8888888888888888d0 + (656*dlnx)/27
     &     .d0 - 32*dlxm1 - (16*xm1)/3.d0 -16*dlnx*xm1 + (64*dlxm1*xm1)
     &     /3.d0 + (8*xm1**2)/9.d0 +(80*dlnx*xm1**2)/27.d0 - (32*dlxm1
     &     *xm1**2)/9.d0)*dli2a +(-0.2962962962962963d0 - (8*xm1**2)
     &     /27.d0)*dli3a +(-32.2962962962963d0 + (64*xm1)/3.d0 - (104
     &     *xm1**2)/27.d0)*dli3d

c..   lfh-terms:
      deltaqqa = deltaqqa + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8*dlnx
     &     *xm1)/3.d0 + (8*xm1**2) /9.d0 -(4*dlnx*xm1**2)/9.d0) +lfh*((4
     &     *dlnx)/9.d0 - 4*dlnx**2 + 16 *dlnx*dlxm1 - (44*xm1)/3.d0 -8
     &     *dlnx*xm1 + (8*dlnx**2*xm1)/3.d0 + (128*dlxm1*xm1)/9.d0 -(32
     &     *dlnx*dlxm1*xm1)/3.d0 + (10*xm1**2)/3.d0 + (20*dlnx*xm1**2)/9
     &     .d0 -(4*dlnx**2*xm1**2)/9.d0 - (32*dlxm1*xm1**2)/9.d0 +(16
     &     *dlnx*dlxm1*xm1**2)/9.d0 +(16 - (32*xm1)/3.d0 + (16*xm1 **2)
     &     /9.d0)*dli2a)

c..   there are no lfr-terms

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaqqa = deltaqqa - ((272*dlnx)/27.d0 - (64*dlnx**2)/27.d0 +
     &        (272*xm1)/27.d0 - (176*dlnx*xm1)/27.d0 + (64*dlnx**2*xm1)
     &        /27.d0 + (8*xm1**2)/9.d0 )
c      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function deltaqqf:

      real*8 function deltaqqf(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

c..   pseudo-scalar
      deltaqqf = 0.d0 * xx

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function deltaqqpa:

      real*8 function deltaqqpa(xx,lfh,lfr)
c..   
c..   qq' contribution (different quarks)
c..   
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar:
      deltaqqpa = 12*dlnx + 4*dlnx**3 - (8*dlnx*dlxm1)/9.d0 + 8*dlnx**2
     &     *dlxm1 -16*dlnx*dlxm1**2 - (152*xm1)/9.d0 - 28*dlnx*xm1 -(32
     &     *dlnx**2*xm1)/9.d0 - (8*dlnx**3*xm1)/3.d0 + (88*dlxm1*xm1)/3
     &     .d0 +16 *dlnx*dlxm1*xm1 - (16*dlnx**2*dlxm1*xm1)/3.d0 -(128
     &     *dlxm1**2 *xm1)/9.d0 + (32*dlnx*dlxm1**2*xm1)/3.d0 + (10*xm1
     &     **2)/3.d0 +(58 *dlnx*xm1**2)/9.d0 + (8*dlnx**2*xm1**2)/9.d0
     &     +(4*dlnx**3*xm1**2) /9.d0 - (20*dlxm1*xm1**2)/3.d0 -(40*dlnx
     &     *dlxm1*xm1**2)/9.d0 + (8 *dlnx**2*dlxm1*xm1**2)/9.d0 +(32
     &     *dlxm1**2*xm1**2)/9.d0 - (16*dlnx *dlxm1**2*xm1**2)/9.d0 +16
     &     *dlnx*z2 + (128*xm1*z2)/9.d0 - (32*dlnx *xm1*z2)/3.d0 -(32
     &     *xm1**2*z2)/9.d0 + (16*dlnx*xm1**2*z2)/9.d0 +(-0
     &     .8888888888888888d0 + 24*dlnx - 32*dlxm1 - (16*xm1)/3.d0 -16
     &     *dlnx *xm1 + (64*dlxm1*xm1)/3.d0 + (8*xm1**2)/9.d0 +(8*dlnx
     &     *xm1**2)/3.d0 - (32*dlxm1*xm1**2)/9.d0)*dli2a +(-32 + (64
     &     *xm1)/3.d0 - (32 *xm1**2)/9.d0)*dli3d

c..   lfh-terms:
      deltaqqpa = deltaqqpa + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8*dlnx
     &     *xm1)/3.d0 + (8*xm1**2) /9.d0 -(4*dlnx*xm1**2)/9.d0) +lfh*((4
     &     *dlnx)/9.d0 - 4*dlnx**2 + 16 *dlnx*dlxm1 - (44*xm1)/3.d0 -8
     &     *dlnx*xm1 + (8*dlnx**2*xm1)/3.d0 + (128*dlxm1*xm1)/9.d0 -(32
     &     *dlnx*dlxm1*xm1)/3.d0 + (10*xm1**2)/3.d0 + (20*dlnx*xm1**2)/9
     &     .d0 -(4*dlnx**2*xm1**2)/9.d0 - (32*dlxm1*xm1**2)/9.d0 +(16
     &     *dlnx*dlxm1*xm1**2)/9.d0 +(16 - (32*xm1)/3.d0 + (16*xm1 **2)
     &     /9.d0)*dli2a)
      
c..   there are no lfr-terms

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaqqpa = deltaqqpa - ((80*dlnx)/9.d0 - (16*dlnx**2)/9.d0 +
     &        (80*xm1)/9.d0 - (16*dlnx *xm1)/3.d0+(16*dlnx**2*xm1)/9.d0
     &        + (8*xm1**2)/9.d0 )
c      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function deltaqqpf:

      real*8 function deltaqqpf(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

c..   pseudo-scalar:
      deltaqqpf = 0.d0 * xx

c..   scalar and pseudo-scalar are the same.

c-- checked against checks.m [17/06/09,rh]

      end

c-}}}
c-{{{ function deltaqqba:

      real*8 function deltaqqba(xx,lfh,lfr)
c..
c..   q q-bar contribution
c..
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

      dli2a = dilog(xm1)
      dli2b = dilog(1.d0-xx**2)
      dli3a = trilog(xm1)
      dli3b = trilog(-(xm1/(2.d0 - xm1)))
      dli3c = trilog(xm1/(2.d0 - xm1))
      dli3d = trilog(-(xm1/xx)) 
      dli3e = trilog(1.d0 - xx**2)
      dli3f = trilog(-((1.d0 - xx**2)/xx**2))

c..   pseudo-scalar
      deltaqqba = (64*dlnx)/9.d0 - (136*dlnx**2)/81.d0 + (88*dlnx**3)/27
     &     .d0 +(184*dlnx*dlxm1)/81.d0 + 8*dlnx**2*dlxm1 - 16*dlnx*dlxm1
     &     **2 -(2020*xm1)/81.d0 - (1708*dlnx*xm1)/81.d0 - (8*dlnx**2
     &     *xm1)/27.d0 -(56*dlnx**3*xm1)/27.d0 + (2632*dlxm1*xm1)/81.d0
     &     +(304*dlnx*dlxm1*xm1)/27.d0 - (16*dlnx**2*dlxm1*xm1)/3.d0
     &     -(128*dlxm1**2*xm1)/9.d0 + (32*dlnx*dlxm1**2*xm1)/3.d0 +(860
     &     *xm1**2)/81.d0 - (386*dlnx*xm1**2)/81.d0 -(136*dlnx**2*xm1**2
     &     )/27.d0 + (8*dlnx**3*xm1**2)/27.d0 -(796*dlxm1*xm1**2)/81.d0
     &     + (8*dlnx*dlxm1*xm1**2)/27.d0 +(8*dlnx**2*dlxm1*xm1**2)/9.d0
     &     + (32*dlxm1**2*xm1**2)/9.d0 -(16*dlnx*dlxm1**2*xm1**2)/9.d0 +
     &     (1060*xm1**3)/27.d0 +(512*dlnx*xm1**3)/27.d0 + (352*dlnx**2
     &     *xm1**3)/81.d0 -(512*dlxm1*xm1**3)/27.d0 - (512*dlnx*dlxm1
     &     *xm1**3)/81.d0 +(416*dlxm1**2*xm1**3)/81.d0 + 16*dlnx*z2 +
     &     (128*xm1*z2)/9.d0 -(32*dlnx*xm1*z2)/3.d0 - (32*xm1**2*z2)/9
     &     .d0 + (16*dlnx*xm1**2*z2)/9.d0 -(688*xm1**3*z2)/81.d0 

      deltaqqba = deltaqqba + (7 .604938271604938d0 + (256*dlnx)/9.d0
     &     -32*dlxm1 - (208*xm1)/9 .d0 - (176*dlnx*xm1)/9.d0 + (64*dlxm1
     &     *xm1)/3.d0 + (152*xm1**2)/9.d0 + (32*dlnx*xm1**2)/9.d0 - (32
     &     *dlxm1*xm1**2)/9.d0 - (208*xm1**3)/81.d0)*dli2a + (-2
     &     .6666666666666665d0 - (20*dlnx)/9.d0 + (176*xm1)/27.d0 + (16
     &     *dlnx*xm1)/9.d0 - (152*xm1**2)/27.d0 - (4*dlnx*xm1**2)/9.d0 +
     &     (16*xm1**3)/9.d0)*dli2b + (-7.407407407407407d0 + (160*xm1)
     &     /27.d0 - (40*xm1**2)/27.d0)*dli3a + (-1 .4814814814814814d0 +
     &     (32*xm1)/27.d0 - (8*xm1**2)/27.d0) * dli3b + (1
     &     .4814814814814814d0 - (32*xm1)/27 .d0 + (8*xm1**2)/27.d0)
     &     * dli3c + (-36 .44444444444444d0 + (224*xm1)/9.d0 - (40*xm1
     &     **2)/9.d0)*dli3d + (1.8518518518518519d0 - (40*xm1)/27.d0 +
     &     (10*xm1 **2)/27.d0)* dli3e + (1.1111111111111112d0 - (8*xm1)
     &     /9.d0 + (2*xm1**2)/9.d0)* dli3f

c..   lfh-terms:
      deltaqqba = deltaqqba + lfh**2*(-4*dlnx - (32*xm1)/9.d0 + (8*dlnx
     &     *xm1)/3.d0 + (8*xm1**2)/9.d0 -(4*dlnx*xm1**2)/9.d0) +lfh*((
     &     -220 *dlnx)/81.d0 - 4*dlnx**2 + 16*dlnx*dlxm1 - (1444*xm1)/81
     &     .d0 -(88 *dlnx*xm1)/27.d0 + (8*dlnx**2*xm1)/3.d0 + (128*dlxm1
     &     *xm1)/9.d0 -(32 *dlnx*dlxm1*xm1)/3.d0 + (526*xm1**2)/81.d0
     &     -(68*dlnx*xm1**2)/27.d0 - (4*dlnx**2*xm1**2)/9.d0 -(32*dlxm1
     &     *xm1**2)/9.d0 + (16*dlnx *dlxm1*xm1**2)/9.d0 +(88*xm1**3)/9
     &     .d0 + (256*dlnx*xm1**3)/81.d0 - (256*dlxm1*xm1**3)/81.d0 +(16
     &     - (32*xm1)/3.d0 + (16*xm1**2)/9.d0) *dli2a)

c..   lfr-terms:
      deltaqqba = deltaqqba + (-88*lfr*xm1**3)/9.d0

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaqqba = deltaqqba - ((352*dlnx)/27.d0 + (32*dlnx**2)/27.d0
     &        +(352*xm1)/27.d0 - (256*dlnx*xm1)/27.d0 - (32*dlnx**2*xm1)
     &        /27.d0- (64*xm1**2)/9.d0 + (16*xm1**3)/27.d0 )
c..   no additional lfr- and lfh-terms
c      endif

c-- checked against checks.m [17/06/09,rh]

      end

C-}}}
c-{{{ function deltaqqbf:

      real*8 function deltaqqbf(xx,lfh,lfr)
      implicit real*8 (a-h,o-z)
      complex*16 dilog,trilog
      external dilog,trilog
      double precision lfh, lfr, xx
      real*8 z2,z3,xm1,dlnx,dlxm1,dli2a,dli2b,dli3a,dli3b,dli3c,dli3d,
     -     dli3e,dli3f

      z2 = 1.6449340668482264366d0
      z3 = 1.2020569031595942854d0

      xm1 = 1.d0 - xx
      dlnx = dlog(xx)
      dlxm1 = dlog(xm1)

c..   pseudo-scalar:
      deltaqqbf = (-32*dlnx)/81. - (32*xm1)/81. + (16*dlnx*xm1)/27. +
     &     (32*xm1**2)/81. - (328*xm1**3)/243. - (64*dlnx*xm1**3)/81. +
     &     (32*dlxm1*xm1**3)/81.

c..   lfh-terms:
      deltaqqbf = deltaqqbf + (-16*lfh*xm1**3)/27.d0

c..   lfr-terms:
      deltaqqbf = deltaqqbf + (16*lfr*xm1**3)/27.d0

c..   turn to scalar:
c      if (.not.lpseudo) then
         deltaqqbf = deltaqqbf - ((-32*dlnx)/27.d0 - (32*xm1)/27.d0 +
     &        (32*dlnx*xm1)/27.d0 + (16 *xm1**2)/27.d0 )
c..   no additional lfh- and lfr-terms
c      endif

c-- checked against checks.m [17/06/09,rh]

      end






      FUNCTION LI2(x)                                                     
      implicit none                                                           
*     !! Dilogarithm for arguments x < = 1.0                             
      real*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO  
      real*8 C(0:18),H,ALFA,B0,B1,B2,LI2OLD                 
      real*8 Li2                                             
      integer  i                                                       
      
      DATA ZERO /0.0d0/, ONE /1.0d0/                               
      DATA HALF /0.5d0/, MALF /-0.5d0/                             
      DATA MONE /-1.0d0/, MTWO /-2.0d0/                            
      DATA PI3 /3.289868133696453d0/, PI6 /1.644934066848226d0/    
                                                                           
      DATA C( 0) / 0.4299669356081370d0/                              
      DATA C( 1) / 0.4097598753307711d0/                              
      DATA C( 2) /-0.0185884366501460d0/                              
      DATA C( 3) / 0.0014575108406227d0/                              
      DATA C( 4) /-0.0001430418444234d0/                              
      DATA C( 5) / 0.0000158841554188d0/                              
      DATA C( 6) /-0.0000019078495939d0/                              
      DATA C( 7) / 0.0000002419518085d0/                              
      DATA C( 8) /-0.0000000319334127d0/                              
      DATA C( 9) / 0.0000000043454506d0/                              
      DATA C(10) /-0.0000000006057848d0/                              
      DATA C(11) / 0.0000000000861210d0/                              
      DATA C(12) /-0.0000000000124433d0/                              
      DATA C(13) / 0.0000000000018226d0/                              
      DATA C(14) /-0.0000000000002701d0/                              
      DATA C(15) / 0.0000000000000404d0/                              
      DATA C(16) /-0.0000000000000061d0/                              
      DATA C(17) / 0.0000000000000009d0/                              
      DATA C(18) /-0.0000000000000001d0/                              
      
      if(x .gt. 1.00000000001d0) then                                    
         write(6,*)'problems in LI2'
         write(6,*)'x=',x 
         stop                                               
      elseif(x .gt. 1.0d0) then                                          
         x = 1.d0                                                      
      endif                                                              
      IF(X .EQ. ONE) THEN                                                
         LI2OLD=PI6
         LI2=LI2OLD                                                       
         RETURN                                                            
      ELSE IF(X .EQ. MONE) THEN                                          
         LI2OLD=MALF*PI6
         LI2=LI2OLD                                                  
         RETURN                                                            
      END IF                                                             
      T=-X                                                               
      IF(T .LE. MTWO) THEN                                               
         Y=MONE/(ONE+T)                                                    
         S=ONE                                                             
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)                        
      ELSE IF(T .LT. MONE) THEN                                          
         Y=MONE-T                                                          
         S=MONE                                                            
         A=LOG(-T)                                                         
         A=-PI6+A*(A+LOG(ONE+ONE/T))                                       
      ELSE IF(T .LE. MALF) THEN                                          
         Y=(MONE-T)/T                                                      
         S=ONE                                                             
         A=LOG(-T)                                                         
         A=-PI6+A*(MALF*A+LOG(ONE+T))                                      
      ELSE IF(T .LT. ZERO) THEN                                          
         Y=-T/(ONE+T)                                                      
         S=MONE                                                            
         A=HALF*LOG(ONE+T)**2                                              
      ELSE IF(T .LE. ONE) THEN                                           
         Y=T                                                               
         S=ONE                                                             
         A=ZERO                                                            
      ELSE                                                               
         Y=ONE/T                                                           
         S=MONE                                                            
         A=PI6+HALF*LOG(T)**2                                              
      END IF                                                             
      
      H=Y+Y-ONE                                                          
      ALFA=H+H                                                           
      B1=ZERO                                                            
      B2=ZERO                                                            
      DO  I = 18,0,-1                                                    
         B0=C(I)+ALFA*B1-B2                                               
         B2=B1                                                            
         B1=B0                                                            
      ENDDO                                                              
      LI2OLD=-(S*(B0-H*B2)+A) 
      LI2=LI2OLD
      end     
c
c                                   
C......Function Li3(x) for -1 < x < 1 (From Daniel)

      function LI3(x)
      implicit none
      double precision LI3,xlog,x,PI,Z3
      PI=3.14159265358979312D0
      Z3=1.20205690315959429D0
      
      if (x.lt.0.5d0) then
         xlog=dlog(1d0-x)
         LI3=  -xlog -(3*xlog**2)/8.-(17*xlog**3)/216.-(5*xlog**4)/576
     -        - (7*xlog**5)/54000. + (7*xlog**6)/86400. + 19*xlog**7/5556600
     -        - xlog**8/752640-11*xlog**9/127008000+11*xlog**10/435456000
      elseif (x.lt.1d0) then
         xlog=dlog(x)
         LI3=Z3+(Pi**2*xlog)/6+(3d0/4-dlog(-xlog)/2)*xlog**2 
     -        -xlog**3/12-xlog**4/288+xlog**6/86400-xlog**8/10160640  
      elseif (x.eq.1d0) then
         LI3=1.20205690315959429D0
      else
         write(6,*)'wrong argument of Li3!!' 
      endif
      return
      end   
          

      FUNCTION S12(z)                                                     
      implicit none
      external LI2,LI3
      real *8 S12,z,Z3,LI2,LI3,t
      Z3=1.20205690315959429D0
c
      if(z.gt.0d0) then
         S12=(dlog(1-z)**2*dlog(z)+2*dlog(1-z)*LI2(1-z)
     -        -2*LI3(1-z)+2*Z3)/2d0
      elseif(z.eq.0d0) then
         S12=0d0
      elseif(z.lt.0d0) then
         t=-z
         S12=dlog(t)*dlog(1+t)**2/2-dlog(1 + t)**3/3- 
     -        dlog(1+t)*LI2(1/(1+t))-LI3(1/(1 + t))+Z3
      endif
      RETURN
      END













      complex*16 function dilog(xx)
c..
c..   Dilogarithm: dilog(x) = Li_2(x)
c..
c     MARCO change
c      complex*16 xx,wgplg
      complex*16 wgplg
      double precision xx

      dilog = wgplg(1,1,xx)
      end

c..   ------------------------------------------------------------

C-}}}
C-{{{ function trilog:

c..   ------------------------------------------------------------

      complex*16 function trilog(xx)
c..
c..   Trilogarithm: trilog(x) = Li_3(x)
c..
c     MARCO change
c      complex*16 xx,wgplg
      complex*16 wgplg
      double precision xx

      trilog = wgplg(2,1,xx)
      end

