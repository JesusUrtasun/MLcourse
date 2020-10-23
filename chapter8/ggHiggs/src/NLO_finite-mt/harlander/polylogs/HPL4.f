C---- Evaluation of rank 4 HPLs with complex argument 

c --- double version (MarcoB)
      double complex function dcli4(x)
      implicit none
      double precision x
      double complex ris, cli4
      ris = cli4(dcmplx(x,0d0))
      dcli4 = ris
      return
      end

      
c --- Real part of HPL4     
      double precision function HPL4real(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4real = dreal(HPL4(n1,n2,n3,n4,x))
      return
      end

c --- Imaginary part of HPL4     
      double precision function HPL4im(n1,n2,n3,n4,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3,n4
      double complex x,HPL4
      x=dcmplx(xr,xi)
      HPL4im = dimag(HPL4(n1,n2,n3,n4,x))
      return
      end
      
      
C---  maps to 18 basis functions 
      
      double complex function HPL4(n1, n2, n3, n4, x)     
      implicit  none 
      double precision pi, zeta2, zeta3, zeta4
      double complex x, ris, myi
c      double complex xcli2,cli3 ! ATM not needed because of the mappings to basis functions
      double complex cli4, dcli4 ! still needed due to cli4(0.5d0)
      double complex basis1,basis2,basis3,basis4,basis5,basis6,basis7
      double complex basis8,basis9,basis10,basis11,basis12,basis13
      double complex basis14,basis15,basis16,basis17,basis18
      double complex basis3_1,basis3_2,basis3_3,basis3_4,basis3_5
      double complex basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
      integer n1,n2,n3,n4,j 

c     #####################################################
c     basis1(x) = cli4(x) 
c     basis2(x) = cli4(-x)
c     basis3(x) = cli4(1-x)
c     basis4(x) = cli4(1/(1+x)) 
c     basis5(x) = cli4(x/(x-1))
c     basis6(x) = cli4(x/(x+1)) 
c     basis7(x) = cli4((1+x)/2) 
c     basis8(x) = cli4((1-x)/2)
c     basis9(x) = cli4((1-x)/(1+x))
c     basis10(x) = cli4((x-1)/(x+1))
c     basis11(x) = cli4(2x/(1+x))
c     basis12(x) = cli4(2x/(x-1)) 
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
c     basis14(x) = cli4(x^2/(x^2-1)) 
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
c     basis16(x) = ch2m2(x) 
c     basis17(x) = ch21m1(x) 
c     basis18(x) = ch21m1(-x) 
c     #####################################################
c     basis3_1(z) = cli3(z) 
c     basis3_2(z) = cli3(-z)
c     basis3_3(z) = cli3(1-z)
c     basis3_4(z) = cli3(1/(1+z)) 
c     basis3_5(z) = cli3((1+z)/2) 
c     basis3_6(z) = cli3((1-z)/2) 
c     basis3_7(z) = cli3((1-z)/(1+z)) 
c     basis3_8(z) = cli3(2z/(z-1))
c     basis2_1(z) = xcli2(z)
c     basis2_2(z) = xcli2(-z)) 
c     basis2_3(z) = xcli2((1-z)/2)
c     #####################################################
c     #####################################################
    
      pi=3.1415926535897932385d0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      
      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      
      ris=dcmplx(10d70,10d70)
      
      if (x.eq.dcmplx(0d0,0d0)) then
         if (j.ne.41) then
            ris = 0d0
         elseif (j.eq.41) then
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL4: "
            print*, "HPL4(",n1,",",n2,",",n3,",",n4,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      else if (x.eq.dcmplx(1d0,0d0)) then
         if (j.le.54.or.j.eq.68) then
            select case (j)
            case(1)             !-1-1-1-1
               ris = dlog(2d0)**4/24d0
            case(2)             !-1-1-10
               ris = -pi**4/90d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/24d0 + dcli4(0.5d0) + dlog(2d0)*zeta3
            case(3)             !-1-1-11
               ris = pi**4/90d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/12d0 - dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/8d0
            case(4)             !-1-10-1
               ris = pi**4/30d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/8d0 - 3*dcli4(0.5d0) 
     &              - (23*dlog(2d0)*zeta3)/8d0
            case(5)             !-1-100
               ris = pi**4/48d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - dlog(2d0)*zeta3
            case(6)             !-1-101
               ris = pi**4/480d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              - (5*dlog(2d0)*zeta3)/8d0
            case(7)             !-1-11-1
               ris = -pi**4/30d0 - (pi**2*dlog(2d0)**2)/8d0 
     &              + dlog(2d0)**4/12d0 + 3*dcli4(0.5d0) 
     &              + (11*dlog(2d0)*zeta3)/4d0
            case(8)             !-1-110
               ris = (-7*pi**4)/288d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0) 
     &              + (13*dlog(2d0)*zeta3)/8d0
            case(9)             !-1-111
               ris = pi**4/720d0 + dlog(2d0)**4/24d0 
     &              -(dlog(2d0)*zeta3)/8d0
            case(10)            !-10-1-1
               ris = -pi**4/30d0 - (pi**2*dlog(2d0)**2)/8d0 
     &              + dlog(2d0)**4/8d0 + 3*dcli4(0.5d0) 
     &              + (11*dlog(2d0)*zeta3)/4d0
            case(11)            !-10-10
               ris = (-11*pi**4)/288d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + 2*dlog(2d0)*zeta3
            case(12)            !-10-11
               ris = (5*pi**4)/144d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) -dlog(2d0)*zeta3
            case(13)            !-100-1
               ris = -pi**4/288d0 + (3*dlog(2d0)*zeta3)/4d0
            case(14)            !-1000
               ris = (-7*pi**4)/720d0
            case(15)            !-1001
               ris = pi**4/60d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - (3*dlog(2d0)*zeta3)/4d0
            case(16)            !-101-1
               ris = (-7*pi**4)/180d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + (13*dlog(2d0)*zeta3)/8d0
            case(17)            !-1010
               ris = (-17*pi**4)/480d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + (3*dlog(2d0)*zeta3)/2d0
            case(18)            !-1011
               ris = pi**4/288d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/24d0 - dcli4(0.5d0) 
     &              + (dlog(2d0)*zeta3)/8d0
            case(19)            !-11-1-1
               ris = pi**4/30d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 3*dcli4(0.5d0) 
     &              - (23*dlog(2d0)*zeta3)/8d0
            case(20)            !-11-10
               ris = pi**4/360d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)*zeta3
            case(21)            !-11-11
               ris = pi**4/1440d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + (dlog(2d0)*zeta3)/4d0
            case(22)            !-110-1
               ris = (11*pi**4)/240d0 + (pi**2*dlog(2d0)**2)/8d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (13*dlog(2d0)*zeta3)/4d0
            case(23)            !-1100
               ris = (19*pi**4)/1440d0 - (3*dlog(2d0)*zeta3)/4d0
            case(24)            !-1101
               ris = (19*pi**4)/1440d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/24d0 - dcli4(0.5d0) 
     &              - (dlog(2d0)*zeta3)/4d0
            case(25)            !-111-1
               ris = -pi**4/288d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + (7*dlog(2d0)*zeta3)/8d0
            case(26)            !-1110
               ris = -pi**4/720d0 - dlog(2d0)**4/24d0 - dcli4(0.5d0) 
     &              + (dlog(2d0)*zeta3)/8d0
            case(27)            !-1111
               ris = dcli4(0.5d0)
            case(28)            !0-1-1-1
               ris = pi**4/90d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/24d0 - dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/8d0
            case(29)            !0-1-10
               ris = -pi**4/288d0
            case(30)            !0-1-11
               ris = -pi**4/80d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0)
            case(31)            !0-10-1
               ris = (13*pi**4)/288d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/2d0
            case(32)            !0-100
               ris = (7*pi**4)/240d0
            case(33)            !0-101
               ris = pi**4/480d0
            case(34)            !0-11-1
               ris = (-7*pi**4)/720d0 - (pi**2*dlog(2d0)**2)/4d0 
     &              + (21*dlog(2d0)*zeta3)/8d0
            case(35)            !0-110
               ris = (13*pi**4)/1440d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0)
            case(36)            !0-111
               ris = (-11*pi**4)/720d0 + dlog(2d0)**4/8d0 
     &              + 3*dcli4(0.5d0)
            case(37)            !00-1-1
               ris = -pi**4/48d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0) 
     &              + (7*dlog(2d0)*zeta3)/4d0
            case(38)            !00-10
               ris = (-7*pi**4)/240d0
            case(39)            !00-11
               ris = -pi**4/180d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0)
            case(40)            !000-1
               ris = (7*pi**4)/720d0
            case(41)            !0000
               ris = 0d0
            case(42)            !0001
               ris = pi**4/90d0
            case(43)            !001-1
               ris = (-19*pi**4)/1440d0 + (7*dlog(2d0)*zeta3)/4d0
            case(44)            !0010
               ris = -pi**4/30d0
            case(45)            !0011
               ris = pi**4/360d0
            case(46)            !01-1-1
               ris = (7*pi**4)/288d0 + (5*pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - (21*dlog(2d0)*zeta3)/8d0
            case(47)            !01-10
               ris = (-11*pi**4)/480d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0)
            case(48)            !01-11
               ris = (7*pi**4)/288d0 - (pi**2*dlog(2d0)**2)/8d0 
     &              - dlog(2d0)**4/8d0 - 3*dcli4(0.5d0)
            case(49)            !010-1
               ris = (71*pi**4)/1440d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/2d0
            case(50)            !0100
               ris = pi**4/30d0
            case(51)            !0101
               ris = pi**4/120d0
            case(52)            !011-1
               ris = -pi**4/80d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/24d0 + dcli4(0.5d0) 
     &              + (7*dlog(2d0)*zeta3)/8d0
            case(53)            !0110
               ris = -pi**4/72d0
            case(54)            !0111
               ris = pi**4/90d0
            case(68)            !1000
               ris = -pi**4/90d0
            end select
         else
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL4: "
            print*, "HPL4(",n1,",",n2,",",n3,",",n4,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      else if (x.eq.dcmplx(-1d0,0d0)) then
         if (j.ge.28) then
            select case (j)
            case(28)            !0-1-1-1
               ris = -pi**4/90d0
            case(29)            !0-1-10
               ris = -pi**4/72d0 + myi*pi*zeta3
            case(30)            !0-1-11
               ris = pi**4/80d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/24d0 - dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/8d0
            case(31)            !0-10-1
               ris = pi**4/120d0
            case(32)            !0-100
               ris = pi**4/20d0 + 2d0*myi*pi*zeta3
            case(33)            !0-101
               ris = (71*pi**4)/1440d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/2d0
            case(34)            !0-11-1
               ris = (-7*pi**4)/288d0 + (pi**2*dlog(2d0)**2)/8d0 
     &              + dlog(2d0)**4/8d0 + 3*dcli4(0.5d0)
            case(35)            !0-110
               ris = (-71*pi**4)/1440d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + myi*pi*((pi**2*dlog(2d0))/4d0 
     &              - zeta3) + (7*dlog(2d0)*zeta3)/2d0 
     &              - 2*((-19*pi**4)/1440d0 
     &              + (7*dlog(2d0)*zeta3)/4d0)
            case(36)            !0-111
               ris = (-7*pi**4)/288d0 - (5*pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0) 
     &              + (21*dlog(2d0)*zeta3)/8d0
            case(37)            !00-1-1
               ris = pi**4/360d0
            case(38)            !00-10
               ris = pi**4/30d0 - myi*pi*zeta3
            case(39)            !00-11
               ris = (-19*pi**4)/1440d0 + (7*dlog(2d0)*zeta3)/4d0
            case(40)            !000-1
               ris = -pi**4/90d0
            case(41)            !0000
               ris = pi**4/24d0
            case(42)            !0001
               ris = (-7*pi**4)/720d0
            case(43)            !001-1
               ris = -pi**4/180d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0)
            case(44)            !0010
               ris = (7*pi**4)/240d0 - myi*3d0/4d0*pi*zeta3
            case(45)            !0011
               ris = -pi**4/48d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0) 
     &              + (7*dlog(2d0)*zeta3)/4d0
            case(46)            !01-1-1
               ris = (11*pi**4)/720d0 - dlog(2d0)**4/8d0 -3*dcli4(0.5d0)
            case(47)            !01-10
               ris = -pi**4/480d0 - 2*(-pi**4/180d0 
     &              - (pi**2*dlog(2d0)**2)/12d0 + dlog(2d0)**4/12d0 
     &              + 2*dcli4(0.5d0)) 
     &              + myi*pi*(-(pi**2*dlog(2d0))/4d0 + (13*zeta3)/8d0)
            case(48)            !01-11
               ris = (7*pi**4)/720d0 + (pi**2*dlog(2d0)**2)/4d0 
     &              - (21*dlog(2d0)*zeta3)/8d0
            case(49)            !010-1
               ris = pi**4/480d0
            case(50)            !0100
               ris = pi**4/80d0 + myi*3d0/2d0*pi*zeta3
            case(51)            !0101
               ris = (13*pi**4)/288d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/2d0
            case(52)            !011-1
               ris = pi**4/80d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0)
            case(53)            !0110
               ris = (-13*pi**4)/288d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + myi/8d0*pi*zeta3 
     &              + (7*dlog(2d0)*zeta3)/2d0 - 2*(-pi**4/48d0 
     &              - (pi**2*dlog(2d0)**2)
     &              /12d0 + dlog(2d0)**4/12d0 + 2*dcli4(0.5d0) 
     &              + (7*dlog(2d0)*zeta3)/4d0)
            case(54)            !0111
               ris = -pi**4/90d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + dcli4(0.5d0) 
     &              + (7*dlog(2d0)*zeta3)/8d0
            case(55)            !1-1-1-1
               ris = dcli4(0.5d0)
            case(56)            !1-1-10
               ris = pi**4/720d0 + dlog(2d0)**4/24d0 + dcli4(0.5d0) 
     &              - myi*pi*(-(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &              + (7*zeta3)/8d0) - (dlog(2d0)*zeta3)/8d0
            case(57)            !1-1-11
               ris = -pi**4/288d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + (7*dlog(2d0)*zeta3)/8d0
            case(58)            !1-10-1
               ris = (-19*pi**4)/1440d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + dcli4(0.5d0) 
     &              + (dlog(2d0)*zeta3)/4d0
            case(59)            !1-100
               ris = (19*pi**4)/1440d0 - (pi**2*(pi**2/12d0 
     &              - dlog(2d0)**2/2d0))/2d0 
     &              - myi*pi*((pi**2*dlog(2d0))/6d0 
     &              - (5*zeta3)/8d0) - (3*dlog(2d0)*zeta3)/4d0 
     &              - myi*pi*(-(pi**2*dlog(2d0))/4d0 + (13*zeta3)/8d0)
            case(60)            !1-101
               ris = (-11*pi**4)/240d0 - (pi**2*dlog(2d0)**2)/8d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + (13*dlog(2d0)*zeta3)/4d0
            case(61)            !1-11-1
               ris = pi**4/1440d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + (dlog(2d0)*zeta3)/4d0
            case(62)            !1-110
               ris = -pi**4/360d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              - myi*pi*((pi**2*dlog(2d0))/12d0 - dlog(2d0)**3/6d0 
     &              - zeta3/4d0) 
     &              + dlog(2d0)*zeta3
            case(63)            !1-111
               ris = pi**4/30d0 + (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/6d0 - 3*dcli4(0.5d0) 
     &              - (23*dlog(2d0)*zeta3)/8d0
            case(64)            !10-1-1
               ris = -pi**4/288d0 - (pi**2*dlog(2d0)**2)/24d0 
     &              + dlog(2d0)**4/24d0 + dcli4(0.5d0) 
     &              - (dlog(2d0)*zeta3)/8d0
            case(65)            !10-10
               ris = -pi**4/480d0 + myi*pi*((pi**2*dlog(2d0))/6d0 
     &              - (5*zeta3)/8d0) - 2*(pi**4/60d0 
     &              + (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - (3*dlog(2d0)*zeta3)/4d0)
            case(66)            !10-11
               ris = (7*pi**4)/180d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/6d0 - 4*dcli4(0.5d0) 
     &              - (13*dlog(2d0)*zeta3)/8d0
            case(67)            !100-1
               ris = pi**4/60d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - (3*dlog(2d0)*zeta3)/4d0
            case(68)            !1000
               ris = (-23*pi**4)/720d0 + myi/6d0*pi**3*dlog(2d0) 
     &              - myi*3d0/4d0*pi*zeta3
            case(69)            !1001
               ris = -pi**4/288d0 + (3*dlog(2d0)*zeta3)/4d0
            case(70)            !101-1
               ris = (-5*pi**4)/144d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) +dlog(2d0)*zeta3
            case(71)            !1010
               ris = (-13*pi**4)/288d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/6d0 + 4*dcli4(0.5d0) 
     &              + myi*pi*((pi**2*dlog(2d0))
     &              /12d0 - zeta3/4d0) + (7*dlog(2d0)*zeta3)/2d0 
     &              - 2*(-pi**4/288d0 
     &              + (3*dlog(2d0)*zeta3)/4d0)
            case(72)            !1011
               ris = pi**4/30d0 + (pi**2*dlog(2d0)**2)/8d0 
     &              - dlog(2d0)**4/8d0 - 3*dcli4(0.5d0) 
     &              - (11*dlog(2d0)*zeta3)/4d0
            case(73)            !11-1-1
               ris = pi**4/720d0 + dlog(2d0)**4/24d0 
     &              - (dlog(2d0)*zeta3)/8d0
            case(74)            !11-10
               ris = (7*pi**4)/288d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - myi*pi*(-dlog(2d0)**3/6d0 
     &              + zeta3/8d0) - (13*dlog(2d0)*zeta3)/8d0
            case(75)            !11-11
               ris = -pi**4/30d0 - (pi**2*dlog(2d0)**2)/8d0 
     &              + dlog(2d0)**4/12d0 + 3*dcli4(0.5d0) 
     &              + (11*dlog(2d0)*zeta3)/4d0
            case(76)            !110-1
               ris = -pi**4/480d0 - (pi**2*dlog(2d0)**2)/12d0 
     &              + (5*dlog(2d0)*zeta3)/8d0
            case(77)            !1100
               ris = pi**4/48d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              - dlog(2d0)**4/12d0 - 2*dcli4(0.5d0) 
     &              - myi*pi*((pi**2*dlog(2d0))
     &              /12d0 - zeta3/4d0) - myi/8d0*pi*zeta3 
     &              - dlog(2d0)*zeta3
            case(78)            !1101
               ris = -pi**4/30d0 - (pi**2*dlog(2d0)**2)/6d0 
     &              + dlog(2d0)**4/8d0 + 3*dcli4(0.5d0) 
     &              + (23*dlog(2d0)*zeta3)/8d0
            case(79)            !111-1
               ris = pi**4/90d0 + (pi**2*dlog(2d0)**2)/24d0 
     &              - dlog(2d0)**4/12d0 - dcli4(0.5d0) 
     &              - (7*dlog(2d0)*zeta3)/8d0
            case(80)            !1110
               ris = pi**4/90d0 + (pi**2*dlog(2d0)**2)/12d0 
     &              - myi/6d0*pi*dlog(2d0)**3 - dlog(2d0)**4/24d0 
     &              - dcli4(0.5d0) 
     &              - dlog(2d0)*zeta3
            case(81)            !1111
               ris = dlog(2d0)**4/24d0
            end select
         else
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL4: "
            print*, "HPL4(",n1,",",n2,",",n3,",",n4,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      else
c ##########  hack to get branch cuts right #################
         if (dimag(x).eq.0d0) then
            x = x + dcmplx(0d0,1d-60)
         endif
c ###########################################################
      select case (j)
            
      case(1)                !-1-1-1-1
            
            ris = log(1d0+x)**4/24d0
            
      case(2)                   !-1-1-10
         
         ris = -pi**4/90d0 + basis4(x) - (pi**2*log(1d0+x)**2)
     &/12d0 + log(1d0+x)**4/24d0 + log(1d0+x)*zeta3
         
      case(3)                   !-1-1-11
         
         ris = -cli4(dcmplx(0.5d0,0d0)) + basis7(x) 
     &+ (pi**2*dlog(2d0)
     &*log(1d0+x))/12d0 - (dlog(2d0)**3*log(1d0+x))/6d0 - (pi**2
     &*log(1d0+x)**2)/24d0 + (dlog(2d0)**2*log(1d0+x)**2)/4d0 
     &- (dlog(2d0)*log(1d0+x)**3)/6d0 - (7*log(1d0+x)*zeta3)/8d0

      case(4)                   !-1-10-1

         ris = pi**4/30d0 - 3*basis4(x) - basis3_4(x)
     &*log(1d0+x) + (pi**2*log(1d0+x)**2)/12d0 + log(1d0+x)**4/24d0 
     &- 2*log(1d0+x)*zeta3

      case(5)                   !-1-100

         ris = pi**4/90d0 - basis2(x) - basis4(x) 
     &- basis6(x) - basis3_4(x)*log(x) - (pi**2*log(x)
     &*log(1d0+x))/6d0 + (pi**2*log(1d0+x)**2)/12d0 - (log(x)**2
     &*log(1d0+x)**2)/4d0 + (log(x)*log(1d0+x)**3)/3d0 
     &- log(1d0+x)**4/12d0 + log(x)*zeta3 - log(1d0+x)*zeta3

      case(6)                   !-1-101

         ris = pi**4/480d0 + basis3(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 - basis13(x)/4d0 
     &+ basis3_4(x)*log(1d0-x) + (pi**2*log(1d0-x)*log(1d0+x))
     &/12d0 + (pi**2*log(1d0+x)**2)/12d0 - (log(1d0-x)*log(1d0+x)**3)
     &/6d0 - (7d0*log(1d0-x)*zeta3)/8d0 - (5d0*log(1d0+x)*zeta3)/8d0
      case(7)                   !-1-11-1

         ris = 3d0*cli4(dcmplx(0.5d0,0d0)) 
     &- 3d0*basis7(x) 
     &+ basis3_5(x)*log(1d0+x) 
     &- (pi**2*dlog(2d0)*log(1d0+x))/6d0 
     &+ (dlog(2d0)**3*log(1d0+x))/3d0 + (pi**2*log(1d0+x)**2)/24d0 
     &- (dlog(2d0)**2*log(1d0+x)**2)/4d0 + (7*log(1d0+x)*zeta3)/4d0

      case(8)                   !-1-110

         ris = pi**4/90d0 + 3*cli4(dcmplx(0.5d0,0d0)) + basis2(x)/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 - basis4(x) 
     &+ basis11(x) - 3*basis7(x) 
     &+ basis3_5(x)*log(x) + (pi**2*dlog(2d0)*log(x))/12d0 
     &- (dlog(2d0)**3*log(x))/6d0 - (pi**2*dlog(2d0)*log(1d0+x))/4d0 
     &+ (dlog(2d0)**3*log(1d0+x))/2d0 - (pi**2*log(x)*log(1d0+x))/12d0 
     &+ (dlog(2d0)**2*log(x)*log(1d0+x))/2d0 + (5*pi**2*log(1d0+x)**2)
     &/24d0 - (3*dlog(2d0)**2*log(1d0+x)**2)/4d0 - (dlog(2d0)*log(x)
     &*log(1d0+x)**2)/2d0 + (dlog(2d0)*log(1d0+x)**3)/2d0 + (log(x)
     &*log(1d0+x)**3)/6d0 - log(1d0+x)**4/6d0 - (7*log(x)*zeta3)/8d0 
     &+ (13*log(1d0+x)*zeta3)/8d0
c     ##############################
         print*, "case 8:"
         print*, basis13(x),basis15(x)
c     ##############################
      case(9)                   !-1-111

         ris = (-7*pi**4)/720d0 - basis8(x) 
     &- basis10(x) + basis7(x) 
     &- basis3_5(x)*log(1d0-x) - (pi**2*dlog(2d0)*log(1d0-x))/6d0 
     &+ (dlog(2d0)**3*log(1d0-x))/3d0 - (dlog(2d0)**2*log(1d0-x)**2)/4d0 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/12d0 -(dlog(2d0)**3*log(1d0+x))/6d0 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/6d0 - (dlog(2d0)**2*log(1d0-x)
     &*log(1d0+x))/2d0 + (dlog(2d0)*log(1d0-x)**2*log(1d0+x))/2d0 
     &- (pi**2*log(1d0+x)**2)/12d0 + (dlog(2d0)**2*log(1d0+x)**2)/4d0 
     &- (log(1d0-x)**2*log(1d0+x)**2)/4d0 + (log(1d0-x)
     &*log(1d0+x)**3)/6d0 - log(1d0+x)**4/24d0 + log(1d0-x)*zeta3 
     &- (log(1d0+x)*zeta3)/8d0

      case(10)                  !-10-1-1

         ris = -pi**4/30d0 + 3*basis4(x) + 2*basis3_4(x)
     &*log(1d0+x) + (pi**2*log(1d0+x)**2)/12d0 + (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (log(x)*log(1d0+x)**3)/2d0 
     &- (5*log(1d0+x)**4)/24d0 + log(1d0+x)*zeta3

      case(11)                  !-10-10

         ris = -pi**4/45d0 + basis2_2(x)**2/2d0 + 2*basis2(x) 
     &+ 2*basis4(x) + 2*basis6(x) + 2*basis3_4(x)
     &*log(x) + (pi**2*log(x)*log(1d0+x))/3d0 + basis2_2(x)*log(x)
     &*log(1d0+x) - (pi**2*log(1d0+x)**2)/6d0 + log(x)**2
     &*log(1d0+x)**2 - (2*log(x)*log(1d0+x)**3)/3d0 + log(1d0+x)**4
     &/6d0 - 2*log(x)*zeta3 + 2*log(1d0+x)*zeta3

      case(12)                  !-10-11

         ris = (-11*pi**4)/720d0 + (pi**2*basis2_2(x))/12d0 
     &- basis2_3(x)*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ basis18(x) 
     &- 2*basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 + basis4(x) 
     &+ basis9(x) - basis10(x) + basis3(x) 
     &+ basis11(x) + basis13(x)/2d0 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &- 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) 
     &+ (basis2_2(x)*dlog(2d0)**2)/2d0
     &+ basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0
     &- (basis2_2(x)*dlog(2d0)**2)/2d0 - 2*basis3_4(x)*log(1d0-x) 
     &+ basis2_2(x)*dlog(2d0)*log(1d0-x) + 2*basis3_1(x)*log(1d0+x) 
     &+ 2*basis3_8(x)*log(1d0+x) + 2*basis3_7(x)
     &*log(1d0+x) - (pi**2*log(1d0-x)*log(1d0+x))/2d0 - basis2_2(x)
     &*log(1d0-x)*log(1d0+x) + dlog(2d0)*log(1d0-x)**2*log(1d0+x) 
     &- (log(1d0-x)**3*log(1d0+x))/3d0 + log(1d0-x)**2*log(x)
     &*log(1d0+x) + (pi**2*log(1d0+x)**2)/8d0 - (basis2_3(x)
     &*log(1d0+x)**2)/2d0 + (basis2_2(x)*log(1d0+x)**2)/2d0
     &-(basis2_1(x)
     &*log(1d0+x)**2)/2d0 - (dlog(2d0)**2*log(1d0+x)**2)/4d0 
     &- (3*dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 - 2*log(1d0-x)
     &*log(x)*log(1d0+x)**2 + dlog(2d0)*log(1d0+x)**3 + (5*log(1d0-x)
     &*log(1d0+x)**3)/6d0 + (5*log(x)*log(1d0+x)**3)/6d0 
     &- (11*log(1d0+x)**4)/24d0 + (7*log(1d0-x)*zeta3)/4d0 
     &+ (log(1d0+x)*zeta3)/4d0

      case(13)                  !-100-1

         ris = -basis2_2(x)**2/2d0 - basis3_2(x)*log(1d0+x)

      case(14)                  !-1000

         ris = basis2(x) - basis3_2(x)*log(x) + (basis2_2(x)
     &*log(x)**2)/2d0 + (log(x)**3*log(1d0+x))/6d0

      case(15)                  !-1001

         ris = -pi**4/360d0 + basis2_2(x)*basis2_1(x) + basis16(x) 
     &- basis3(x) - basis2(x) - basis1(x) + basis5(x) 
     &+ basis4(x) + basis6(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*log(1d0-x) 
     &+ (pi**2*log(1d0-x)**2)/16d0 + log(1d0-x)**4/32d0 
     &- (log(1d0-x)**3*log(x))/12d0 + 2*basis3_1(x)*log(1d0+x) 
     &- (pi**2*log(1d0-x)*log(1d0+x))/24d0 - (log(1d0-x)**3
     &*log(1d0+x))/24d0 + (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &- (5*pi**2*log(1d0+x)**2)/48d0 - (log(1d0-x)**2*log(1d0+x)**2)
     &/16d0 + (log(1d0-x)*log(x)*log(1d0+x)**2)/4d0 - (log(1d0-x)
     &*log(1d0+x)**3)/24d0 - (log(x)*log(1d0+x)**3)/12d0 
     &+ (7*log(1d0+x)**4)/96d0 + (3*log(1d0-x)*zeta3)/4d0 
     &+ (3*log(1d0+x)*zeta3)/4d0

      case(16)                  !-101-1

         ris = pi**4/90d0 - (pi**2*basis2_2(x))/12d0 + basis2_3(x)
     &*basis2_2(x) - basis2_2(x)**2/2d0 
     &- basis18(x) 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0)
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 
     &+ basis15(x)/4d0 - basis4(x) 
     &+ 2*basis6(x) - basis11(x) + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 - basis2_2(x)*dlog(2d0)*log(1d0-x) 
     &+ basis3_6(x)*log(1d0+x) - basis3_3(x)*log(1d0+x) 
     &- 2*basis3_1(x)*log(1d0+x) - 2*basis3_8(x)*log(1d0+x) 
     &- basis3_4(x)*log(1d0+x) - basis3_7(x)
     &*log(1d0+x) + basis3_5(x)*log(1d0+x) + (pi**2*dlog(2d0)
     &*log(1d0+x))/6d0 - (dlog(2d0)**3*log(1d0+x))/3d0 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/4d0 + (dlog(2d0)**2*log(1d0-x)
     &*log(1d0+x))/2d0 - dlog(2d0)*log(1d0-x)**2*log(1d0+x) 
     &+ (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2*log(x)
     &*log(1d0+x) - (3*pi**2*log(1d0+x)**2)/8d0 + (basis2_3(x)
     &*log(1d0+x)**2)/2d0 - (basis2_2(x)*log(1d0+x)**2)/2d0
     &+(basis2_1(x)
     &*log(1d0+x)**2)/2d0 + (3*dlog(2d0)**2*log(1d0+x)**2)/4d0 
     &+ (dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 + log(1d0-x)*log(x)
     &*log(1d0+x)**2 - dlog(2d0)*log(1d0+x)**3 - (5*log(x)
     &*log(1d0+x)**3)/6d0 + (11*log(1d0+x)**4)/24d0 
     &+ (log(1d0+x)*zeta3)/4d0

      case(17)                  !-1010

         ris = -(basis2_2(x)*basis2_1(x)) - basis16(x) 
     &+ basis3_6(x)
     &*log(x) - basis3_3(x)*log(x) - basis3_4(x)*log(x) 
     &+ basis3_7(x)*log(x) + basis3_5(x)*log(x) 
     &+ (pi**2*dlog(2d0)*log(x))/6d0 - (dlog(2d0)**3*log(x))/3d0 
     &- (pi**2*log(1d0-x)*log(x))/12d0 - basis2_2(x)*log(1d0-x)
     &*log(x) + (dlog(2d0)**2*log(1d0-x)*log(x))/2d0 - 2*basis3_1(x)
     &*log(1d0+x) - (pi**2*log(x)*log(1d0+x))/12d0 
     &+ (dlog(2d0)**2*log(x)*log(1d0+x))/2d0 - dlog(2d0)*log(1d0-x)
     &*log(x)*log(1d0+x) - log(1d0-x)*log(x)**2*log(1d0+x) 
     &+ (log(1d0-x)*log(x)*log(1d0+x)**2)/2d0 - (3*log(x)*zeta3)/4d0

      case(18)                  !-1011

         ris = pi**4/288d0 - basis4(x) + basis9(x)
     &/2d0 - basis10(x)/2d0 - basis13(x)/4d0 
     &- basis3_6(x)*log(1d0-x) + basis3_3(x)*log(1d0-x) 
     &+ basis3_4(x)*log(1d0-x) - basis3_7(x)
     &*log(1d0-x) - basis3_5(x)*log(1d0-x) - (pi**2*dlog(2d0)
     &*log(1d0-x))/6d0 + (dlog(2d0)**3*log(1d0-x))/3d0 + (pi**2
     &*log(1d0-x)**2)/24d0 + (basis2_2(x)*log(1d0-x)**2)/2d0 
     &- (dlog(2d0)**2*log(1d0-x)**2)/2d0 + (pi**2*log(1d0-x)
     &*log(1d0+x))/4d0 - (dlog(2d0)**2*log(1d0-x)*log(1d0+x))/2d0 
     &+ dlog(2d0)*log(1d0-x)**2*log(1d0+x) + (log(1d0-x)**2*log(x)
     &*log(1d0+x))/2d0 + (pi**2*log(1d0+x)**2)/24d0 
     &- (log(1d0-x)**2*log(1d0+x)**2)/2d0 - log(1d0+x)**4/24d0 
     &+ (log(1d0-x)*zeta3)/8d0 + (log(1d0+x)*zeta3)/8d0

      case(19)                  !-11-1-1

         ris = -3*cli4(dcmplx(0.5d0,0d0)) + 3*basis7(x) 
     &- 2*basis3_5(x)*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))
     &/12d0 - (dlog(2d0)**3*log(1d0+x))/6d0 + (pi**2*log(1d0+x)**2)/12d0 
     &- (basis2_3(x)*log(1d0+x)**2)/2d0 - (dlog(2d0)**2
     &*log(1d0+x)**2)/2d0 + (dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 
     &+ (dlog(2d0)*log(1d0+x)**3)/2d0 - (log(1d0-x)*log(1d0+x)**3)/2d0 
     &- (7*log(1d0+x)*zeta3)/8d0

      case(20)                  !-11-10

         ris = -pi**4/90d0 - basis2_2(x)**2/2d0 
     &- basis18(x) 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0)
     &- 6*cli4(dcmplx(0.5d0,0d0)) + basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + (3*basis15(x))/4d0 
     &+ basis4(x) + 2*basis6(x) - 3*basis11(x) 
     &+ 6*basis7(x) - 2*basis3_5(x)*log(x) 
     &- (pi**2*dlog(2d0)*log(x))/6d0 + (dlog(2d0)**3*log(x))/3d0 
     &- 2*basis3_1(x)*log(1d0+x) - 2*basis3_8(x)*log(1d0+x) 
     &- 2*basis3_7(x)*log(1d0+x) + (pi**2*dlog(2d0)
     &*log(1d0+x))/2d0 - dlog(2d0)**3*log(1d0+x) 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/3d0 - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 + (pi**2*log(x)
     &*log(1d0+x))/4d0 - basis2_3(x)*log(x)*log(1d0+x) 
     &- (3*dlog(2d0)**2*log(x)*log(1d0+x))/2d0 + dlog(2d0)*log(1d0-x)
     &*log(x)*log(1d0+x) - log(1d0-x)**2*log(x)*log(1d0+x) 
     &- (17*pi**2*log(1d0+x)**2)/24d0 + (basis2_3(x)
     &*log(1d0+x)**2)/2d0 - (basis2_2(x)*log(1d0+x)**2)/2d0 
     &+ (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ (7*dlog(2d0)**2*log(1d0+x)**2)
     &/4d0 + (3*dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 + dlog(2d0)
     &*log(x)*log(1d0+x)**2 + log(1d0-x)*log(x)*log(1d0+x)**2 
     &- 2*dlog(2d0)*log(1d0+x)**3 - (log(1d0-x)*log(1d0+x)**3)/2d0 
     &- (7*log(x)*log(1d0+x)**3)/6d0 + (19*log(1d0+x)**4)/24d0 
     &+ (7*log(x)*zeta3)/4d0 - (9*log(1d0+x)*zeta3)/4d0

      case(21)                  !-11-11

         ris = (11*pi**4)/480d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 + 2*basis8(x) + 2
     &*basis10(x) - 2*basis7(x) 
     &- (pi**2*dlog(2d0)**2)/24d0 + (basis2_3(x)*dlog(2d0)**2)/2d0 
     &+ dlog(2d0)**4/8d0 + 2*basis3_5(x)*log(1d0-x) + (5*pi**2
     &*dlog(2d0)*log(1d0-x))/12d0 - basis2_3(x)*dlog(2d0)
     &*log(1d0-x) - (7*dlog(2d0)**3*log(1d0-x))/6d0 + dlog(2d0)**2
     &*log(1d0-x)**2 - (pi**2*dlog(2d0)*log(1d0+x))/6d0 + (dlog(2d0)**3
     &*log(1d0+x))/3d0 - (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_3(x)*log(1d0-x)*log(1d0+x) + (3*dlog(2d0)**2
     &*log(1d0-x)*log(1d0+x))/2d0 - 2*dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (pi**2*log(1d0+x)**2)/6d0 - (dlog(2d0)**2
     &*log(1d0+x)**2)/2d0 + log(1d0-x)**2*log(1d0+x)**2 
     &- (log(1d0-x)*log(1d0+x)**3)/3d0 + log(1d0+x)**4/12d0 
     &- 2*log(1d0-x)*zeta3 + (log(1d0+x)*zeta3)/4d0

      case(22)                  !-110-1

         ris = -pi**4/90d0 + basis2_2(x)**2/2d0 
     &+ basis18(x) 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0
     &- (3*basis2(x))/2d0 - basis1(x)/2d0 
     &- basis15(x)/4d0 + basis4(x) 
     &- 2*basis6(x) + basis11(x) + basis3_2(x)
     &*log(1d0+x) + basis3_1(x)*log(1d0+x) + basis3_8(x)
     &*log(1d0+x) + basis3_4(x)*log(1d0+x) 
     &+ basis3_7(x)*log(1d0+x) - basis3_5(x)
     &*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/12d0 + (dlog(2d0)**3
     &*log(1d0+x))/6d0 - (pi**2*log(1d0-x)*log(1d0+x))/6d0 
     &+ (dlog(2d0)*log(1d0-x)**2*log(1d0+x))/2d0 - (log(1d0-x)**3
     &*log(1d0+x))/6d0 + (log(1d0-x)**2*log(x)*log(1d0+x))/2d0 
     &+ (3*pi**2*log(1d0+x)**2)/8d0 - (basis2_3(x)
     &*log(1d0+x)**2)/2d0 + (basis2_2(x)*log(1d0+x)**2)/2d0 
     &- (basis2_1(x)*log(1d0+x)**2)/2d0 
     &- (3*dlog(2d0)**2*log(1d0+x)**2)
     &/4d0 - (dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 - log(1d0-x)
     &*log(x)*log(1d0+x)**2 + dlog(2d0)*log(1d0+x)**3 + (5*log(x)
     &*log(1d0+x)**3)/6d0 - (11*log(1d0+x)**4)/24d0 
     &- (log(1d0+x)*zeta3)/8d0

      case(23)                  !-1100

         ris = -pi**4/90d0 - 2*cli4(dcmplx(0.5d0,0d0)) - basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + basis15(x)/4d0 
     &+ basis4(x) + basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis3_2(x)*log(x) - basis3_1(x)*log(x) 
     &- basis3_8(x)*log(x) + basis3_4(x)*log(x) 
     &- basis3_7(x)*log(x) - basis3_5(x)*log(x) 
     &- (pi**2*dlog(2d0)*log(x))/12d0 + (dlog(2d0)**3*log(x))/6d0 
     &+ (pi**2*log(1d0-x)*log(x))/6d0 - (dlog(2d0)*log(1d0-x)**2
     &*log(x))/2d0 + (log(1d0-x)**3*log(x))/6d0 + (pi**2*log(x)**2)
     &/24d0 - (basis2_3(x)*log(x)**2)/2d0 - (dlog(2d0)**2
     &*log(x)**2)/4d0 + (dlog(2d0)*log(1d0-x)*log(x)**2)/2d0 
     &- (log(1d0-x)**2*log(x)**2)/2d0 + (pi**2*dlog(2d0)*log(1d0+x))
     &/6d0 - (dlog(2d0)**3*log(1d0+x))/3d0 + (pi**2*log(x)*log(1d0+x))
     &/12d0 - (dlog(2d0)**2*log(x)*log(1d0+x))/2d0 +dlog(2d0)*log(1d0-x)
     &*log(x)*log(1d0+x) + (log(1d0-x)*log(x)**2*log(1d0+x))/2d0 
     &- (pi**2*log(1d0+x)**2)/6d0 + (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- (log(1d0-x)*log(x)*log(1d0+x)**2)/2d0 - (dlog(2d0)
     &*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/6d0 
     &+ log(1d0+x)**4/6d0 + (7*log(x)*zeta3)/8d0 
     &        - (3*log(1d0+x)*zeta3)/4d0

      case(24)                  !-1101

         ris = pi**4/144d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0)
     &- 4*cli4(dcmplx(0.5d0,0d0)) + 2*basis8(x) - basis3(x) 
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 - 2*basis5(x) 
     &- basis12(x) + basis15(x)/4d0 
     &+ 2*basis4(x) - basis9(x) 
     &+ basis10(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &+ basis14(x)/4d0 - basis3_2(x)
     &*log(1d0-x) + basis3_1(x)*log(1d0-x) - basis3_8(x)
     &*log(1d0-x) - basis3_4(x)*log(1d0-x) 
     &+ basis3_7(x)*log(1d0-x) + basis3_5(x)
     &*log(1d0-x) + (pi**2*dlog(2d0)*log(1d0-x))/4d0 - (dlog(2d0)**3
     &*log(1d0-x))/2d0 - (5*pi**2*log(1d0-x)**2)/48d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (3*dlog(2d0)**2*log(1d0-x)**2)/4d0-(dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (3*log(1d0-x)**4)/32d0 + (log(1d0-x)**3*log(x))/4d0 
     &- 2*basis3_1(x)*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/6d0 
     &- (dlog(2d0)**3*log(1d0+x))/3d0 - (3*pi**2*log(1d0-x)*log(1d0+x))
     &/8d0 + (dlog(2d0)**2*log(1d0-x)*log(1d0+x))/2d0 - dlog(2d0)
     &*log(1d0-x)**2*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/24d0 
     &- (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &- (7*pi**2*log(1d0+x)**2)/48d0 + (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &+ (9*log(1d0-x)**2*log(1d0+x)**2)/16d0 - (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 - (dlog(2d0)*log(1d0+x)**3)/3d0 
     &+ (log(1d0-x)*log(1d0+x)**3)/24d0 - (log(x)*log(1d0+x)**3)
     &/12d0 + (17*log(1d0+x)**4)/96d0 - (log(1d0-x)*zeta3)/8d0 
     &- (7*log(1d0+x)*zeta3)/4d0

      case(25)                  !-111-1

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*dlog(2d0)**2)/24d0 
     &- (basis2_3(x)*dlog(2d0)**2)/2d0 - dlog(2d0)**4/8d0 
     &- (pi**2*dlog(2d0)*log(1d0-x))/12d0 + basis2_3(x)*dlog(2d0)
     &*log(1d0-x) + (dlog(2d0)**3*log(1d0-x))/2d0 - (dlog(2d0)**2
     &*log(1d0-x)**2)/2d0 - basis3_6(x)*log(1d0+x) - (pi**2
     &*dlog(2d0)*log(1d0+x))/12d0 + (dlog(2d0)**3*log(1d0+x))/6d0 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/12d0 - (dlog(2d0)**2*log(1d0-x)
     &*log(1d0+x))/2d0 + (dlog(2d0)*log(1d0-x)**2*log(1d0+x))/2d0 
     &+ (7*log(1d0+x)*zeta3)/8d0

      case(26)                  !-1110

         ris = pi**4/288d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &- (basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0))
     &+ 5*cli4(dcmplx(0.5d0,0d0)) 
     &- basis8(x) + 2*basis3(x) - basis2(x) - basis1(x) 
     &+ 2*basis12(x) - basis15(x)/2d0 
     &- 3*basis4(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 - 2*basis6(x) 
     &+ 4*basis11(x) - 4*basis7(x) 
     &- basis13(x)/4d0 + 2*basis3_8(x)*log(1d0-x) 
     &- (pi**2*dlog(2d0)*log(1d0-x))/12d0 +(dlog(2d0)**3*log(1d0-x))/6d0 
     &- (pi**2*log(1d0-x)**2)/8d0 - (basis2_3(x)
     &*log(1d0-x)**2)/2d0 + (basis2_2(x)*log(1d0-x)**2)/2d0 
     &- (basis2_1(x)*log(1d0-x)**2)/2d0
     &-(dlog(2d0)**2*log(1d0-x)**2)/2d0 
     &+ (2*dlog(2d0)*log(1d0-x)**3)/3d0 - (7*log(1d0-x)**4)/24d0 
     &- basis3_6(x)*log(x) - (pi**2*dlog(2d0)*log(x))/12d0 
     &+ (dlog(2d0)**3*log(x))/6d0 + basis2_3(x)*log(1d0-x)
     &*log(x) - (dlog(2d0)*log(1d0-x)**2*log(x))/2d0 + (log(1d0-x)**3
     &*log(x))/3d0 + 2*basis3_1(x)*log(1d0+x) - (pi**2*dlog(2d0)
     &*log(1d0+x))/3d0 + (2*dlog(2d0)**3*log(1d0+x))/3d0 + (pi**2
     &*log(1d0-x)*log(1d0+x))/6d0 + (3*pi**2*log(1d0+x)**2)/8d0 
     &- dlog(2d0)**2*log(1d0+x)**2 + (2*dlog(2d0)*log(1d0+x)**3)/3d0 
     &+ (log(x)*log(1d0+x)**3)/3d0 - (3*log(1d0+x)**4)/8d0 
     &- (7*log(1d0-x)*zeta3)/4d0 + (7*log(x)*zeta3)/8d0 
     &+ (13*log(1d0+x)*zeta3)/8d0

      case(27)                  !-1111

         ris = cli4(dcmplx(0.5d0,0d0)) - basis8(x) 
     &+ basis3_6(x)
     &*log(1d0-x) - (basis2_3(x)*log(1d0-x)**2)/2d0 + (dlog(2d0)
     &*log(1d0-x)**3)/6d0 - (log(1d0-x)**3*log(1d0+x))/6d0

      case(28)                  !0-1-1-1

         ris = pi**4/90d0 - basis4(x) - basis3_4(x)
     &*log(1d0+x) - (pi**2*log(1d0+x)**2)/12d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ log(1d0+x)**4/8d0

      case(29)                  !0-1-10

         ris = -basis2_2(x)**2/2d0 - basis3_4(x)*log(x) 
     &- (pi**2*log(x)*log(1d0+x))/6d0 - basis2_2(x)*log(x)*log(1d0+x) 
     &- (log(x)**2*log(1d0+x)**2)/2d0 + (log(x)*log(1d0+x)**3)/6d0 
     &+ log(x)*zeta3

      case(30)                  !0-1-11

         ris = -basis18(x)

      case(31)                  !0-10-1

         ris = pi**4/45d0 + basis2_2(x)**2/2d0 - 2*basis2(x) 
     &- 2*basis4(x) - 2*basis6(x) + 2*basis3_2(x)
     &*log(1d0+x) + (pi**2*log(1d0+x)**2)/6d0 + (log(x)
     &*log(1d0+x)**3)/3d0 - log(1d0+x)**4/6d0 - 2*log(1d0+x)*zeta3

      case(32)                  !0-100

         ris = -3*basis2(x) + 2*basis3_2(x)*log(x) 
     &- (basis2_2(x)*log(x)**2)/2d0

      case(33)                  !0-101

         ris = pi**4/180d0 - basis2_2(x)*basis2_1(x) - basis16(x) 
     &+ 2*basis3(x) + 2*basis2(x) + 2*basis1(x) - 2*basis5(x) 
     &- 2*basis4(x) - 2*basis6(x) - basis13(x)/2d0 
     &+ basis14(x)/2d0 - 2*basis3_2(x)*log(1d0-x) 
     &- (pi**2*log(1d0-x)**2)/8d0 - log(1d0-x)**4/16d0 
     &+ (log(1d0-x)**3*log(x))/6d0 - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/12d0 + (log(1d0-x)**3
     &*log(1d0+x))/12d0 - (log(1d0-x)**2*log(x)*log(1d0+x))/2d0 
     &+ (5*pi**2*log(1d0+x)**2)/24d0 + (log(1d0-x)**2*log(1d0+x)**2)
     &/8d0 - (log(1d0-x)*log(x)*log(1d0+x)**2)/2d0 + (log(1d0-x)
     &*log(1d0+x)**3)/12d0 + (log(x)*log(1d0+x)**3)/6d0 
     &- (7*log(1d0+x)**4)/48d0 - (3*log(1d0-x)*zeta3)/2d0 
     &- (3*log(1d0+x)*zeta3)/2d0

      case(34)                  !0-11-1

         ris = pi**4/90d0 + (pi**2*basis2_2(x))/12d0 - basis2_3(x)
     &*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ basis18(x) 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0
     &+ 6*cli4(dcmplx(0.5d0,0d0)) - basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - (3*basis15(x))/4d0 
     &- basis4(x) - 2*basis6(x) + 3*basis11(x) 
     &- 6*basis7(x) - (basis2_2(x)*dlog(2d0)**2)/2d0 + basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - basis3_6(x)*log(1d0+x) 
     &+ basis3_3(x)*log(1d0+x) - basis3_2(x)*log(1d0+x) 
     &+ 3*basis3_1(x)
     &*log(1d0+x) + 3*basis3_8(x)*log(1d0+x) 
     &+ 2*basis3_7(x)*log(1d0+x) - (7*pi**2*dlog(2d0)
     &*log(1d0+x))/12d0 + (7*dlog(2d0)**3*log(1d0+x))/6d0 
     &- (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 - (dlog(2d0)**2*log(1d0-x)
     &*log(1d0+x))/2d0 + (3*dlog(2d0)*log(1d0-x)**2*log(1d0+x))/2d0 
     &- (log(1d0-x)**3*log(1d0+x))/2d0 + (3*log(1d0-x)**2*log(x)
     &*log(1d0+x))/2d0 + (17*pi**2*log(1d0+x)**2)/24d0 
     &- (basis2_3(x)*log(1d0+x)**2)/2d0 + (basis2_2(x)
     &*log(1d0+x)**2)/2d0 - (basis2_1(x)*log(1d0+x)**2)/2d0 
     &- (7*dlog(2d0)**2*log(1d0+x)**2)/4d0 - (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 - 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &+ 2*dlog(2d0)*log(1d0+x)**3 + (log(1d0-x)*log(1d0+x)**3)/2d0 
     &+ (7*log(x)*log(1d0+x)**3)/6d0 - (19*log(1d0+x)**4)/24d0 
     &+ (17*log(1d0+x)*zeta3)/8d0

      case(35)                  !0-110

         ris = pi**4/45d0 + basis2_2(x)*basis2_1(x) + basis16(x) 
     &+ 4*cli4(dcmplx(0.5d0,0d0)) + basis2(x) - 3*basis1(x) 
     &- basis15(x)/2d0 - 2*basis4(x) 
     &- 2*basis6(x) + 4*basis11(x) 
     &- 4*basis7(x) - basis3_6(x)*log(x) 
     &+ basis3_3(x)*log(x) - basis3_2(x)*log(x) 
     &+ basis3_1(x)*log(x) 
     &+ basis3_8(x)*log(x) - (pi**2*dlog(2d0)*log(x))/12d0 
     &+ (dlog(2d0)**3*log(x))/6d0 - (pi**2*log(1d0-x)*log(x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(x) 
     &- (dlog(2d0)**2*log(1d0-x)*log(x))
     &/2d0 + (dlog(2d0)*log(1d0-x)**2*log(x))/2d0 - (log(1d0-x)**3
     &*log(x))/6d0 + (log(1d0-x)**2*log(x)**2)/2d0 + 2*basis3_1(x)
     &*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/3d0 + (2*dlog(2d0)**3
     &*log(1d0+x))/3d0 + (pi**2*log(1d0+x)**2)/3d0 - dlog(2d0)**2
     &*log(1d0+x)**2 + (2*dlog(2d0)*log(1d0+x)**3)/3d0 + (log(x)
     &*log(1d0+x)**3)/3d0 - log(1d0+x)**4/3d0 - (log(x)*zeta3)/8d0 
     &+ (3*log(1d0+x)*zeta3)/2d0

      case(36)                  !0-111

         ris = -pi**4/72d0 - cli4(dcmplx(0.5d0,0d0)) - basis8(x) 
     &- basis3(x) - basis2(x)/2d0 + basis1(x)/2d0 
     &+ 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*log(1d0-x) 
     &- basis3_3(x)*log(1d0-x) + basis3_2(x)*log(1d0-x) 
     &- basis3_1(x)*log(1d0-x) - basis3_8(x)*log(1d0-x) 
     &+ (3*pi**2*log(1d0-x)**2)/16d0 - (basis2_2(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/4d0 - (dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (19*log(1d0-x)**4)/96d0 - (7*log(1d0-x)**3*log(x))/12d0 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/6d0 - (dlog(2d0)**3*log(1d0+x))/3d0 
     &- (pi**2*log(1d0-x)*log(1d0+x))/24d0 - (log(1d0-x)**3
     &*log(1d0+x))/24d0 + (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &- (13*pi**2*log(1d0+x)**2)/48d0 + (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- (log(1d0-x)**2*log(1d0+x)**2)/16d0 + (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 - (dlog(2d0)*log(1d0+x)**3)/3d0 - (log(1d0-x)
     &*log(1d0+x)**3)/24d0 - (log(x)*log(1d0+x)**3)/4d0 
     &+ (23*log(1d0+x)**4)/96d0 + (7*log(1d0-x)*zeta3)/4d0

      case(37)                  !00-1-1

         ris = -pi**4/90d0 + basis2(x) + basis4(x) 
     &+ basis6(x) - basis3_2(x)*log(1d0+x) 
     &- (pi**2*log(1d0+x)**2)/12d0 - (log(x)*log(1d0+x)**3)/6d0 
     &+ log(1d0+x)**4/12d0 + log(1d0+x)*zeta3

      case(38)                  !00-10

         ris = 3*basis2(x) - basis3_2(x)*log(x)

      case(39)                  !00-11

         ris = -pi**4/72d0 - 2*cli4(dcmplx(0.5d0,0d0)) - basis3(x) 
     &- (3*basis2(x))/2d0 + basis1(x)/2d0 + basis5(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*log(1d0-x) 
     &+ (pi**2*log(1d0-x)**2)/16d0 + log(1d0-x)**4/32d0 
     &- (log(1d0-x)**3*log(x))/12d0 + (pi**2*dlog(2d0)*log(1d0+x))/6d0 
     &- (dlog(2d0)**3*log(1d0+x))/3d0 - (pi**2*log(1d0-x)
     &*log(1d0+x))/24d0 - (log(1d0-x)**3*log(1d0+x))/24d0 
     &+ (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 - (13*pi**2
     &*log(1d0+x)**2)/48d0 + (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- (log(1d0-x)**2*log(1d0+x)**2)/16d0 + (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 - (dlog(2d0)*log(1d0+x)**3)/3d0 
     &- (log(1d0-x)*log(1d0+x)**3)/24d0 - (log(x)*log(1d0+x)**3)/4d0 
     &+ (23*log(1d0+x)**4)/96d0 + (3*log(1d0-x)*zeta3)/4d0

      case(40)                  !000-1

         ris = -basis2(x)

      case(41)                  !0000

         ris = log(x)**4/24d0

      case(42)                  !0001

         ris = basis1(x)

      case(43)                  !001-1

         ris = pi**4/90d0 + 2*cli4(dcmplx(0.5d0,0d0)) + basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - basis15(x)/4d0 
     &- basis4(x) - basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) + basis3_1(x)*log(1d0+x) - (pi**2*dlog(2d0)
     &*log(1d0+x))/6d0 + (dlog(2d0)**3*log(1d0+x))/3d0 
     &+ (pi**2*log(1d0+x)**2)/6d0 - (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &+ (dlog(2d0)*log(1d0+x)**3)/3d0 + (log(x)*log(1d0+x)**3)/6d0 
     &- log(1d0+x)**4/6d0 + (3*log(1d0+x)*zeta3)/4d0

      case(44)                  !0010

         ris = -3*basis1(x) + basis3_1(x)*log(x)

      case(45)                  !0011
         
         ris = pi**4/90d0 - basis3(x) + basis1(x) + basis5(x) 
     &- basis3_1(x)*log(1d0-x) + (pi**2*log(1d0-x)**2)/12d0 
     &+ log(1d0-x)**4/24d0 - (log(1d0-x)**3*log(x))/6d0 
     &+ log(1d0-x)*zeta3

      case(46)                  !01-1-1
         ris = -pi**4/90d0 - 3*cli4(dcmplx(0.5d0,0d0)) - basis2(x)/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 + basis4(x) 
     &- basis11(x) + 3*basis7(x) + basis3_2(x)
     &*log(1d0+x) - basis3_1(x)*log(1d0+x) - basis3_8(x)
     &*log(1d0+x) + basis3_4(x)*log(1d0+x) 
     &- basis3_7(x)*log(1d0+x) - basis3_5(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/6d0 - (dlog(2d0)**3
     &*log(1d0+x))/3d0 + (pi**2*log(1d0-x)*log(1d0+x))/6d0 
     &- (dlog(2d0)*log(1d0-x)**2*log(1d0+x))/2d0 + (log(1d0-x)**3
     &*log(1d0+x))/6d0 - (log(1d0-x)**2*log(x)*log(1d0+x))/2d0 
     &- (pi**2*log(1d0+x)**2)/8d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0+x)**2)/4d0 + dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2 + log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (dlog(2d0)*log(1d0+x)**3)/2d0 - (log(1d0-x)*log(1d0+x)**3)/2d0 
     &- (log(x)*log(1d0+x)**3)/6d0 + log(1d0+x)**4/6d0 
     &- (3*log(1d0+x)*zeta3)/4d0

      case(47)                  !01-10
         
         ris = -pi**4/45d0 - basis16(x) - 4*cli4(dcmplx(0.5d0,0d0)) 
     &- basis2(x) 
     &+ 3*basis1(x) + basis15(x)/2d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 4*basis11(x) 
     &+ 4*basis7(x) + basis3_2(x)*log(x) - basis3_1(x)*log(x) 
     &- basis3_8(x)*log(x) + basis3_4(x)*log(x) 
     &- basis3_7(x)*log(x) - basis3_5(x)*log(x) 
     &- (pi**2*dlog(2d0)*log(x))/12d0 + (dlog(2d0)**3*log(x))/6d0 
     &+ (pi**2*log(1d0-x)*log(x))/6d0 - (dlog(2d0)*log(1d0-x)**2
     &*log(x))/2d0 + (log(1d0-x)**3*log(x))/6d0 - (log(1d0-x)**2
     &*log(x)**2)/2d0 - 2*basis3_1(x)*log(1d0+x) + (pi**2*dlog(2d0)
     &*log(1d0+x))/3d0 - (2*dlog(2d0)**3*log(1d0+x))/3d0 + (pi**2
     &*log(x)*log(1d0+x))/12d0 + basis2_1(x)*log(x)*log(1d0+x) 
     &- (dlog(2d0)**2*log(x)*log(1d0+x))/2d0 + dlog(2d0)*log(1d0-x)
     &*log(x)*log(1d0+x) + log(1d0-x)*log(x)**2*log(1d0+x) 
     &- (pi**2*log(1d0+x)**2)/3d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (log(1d0-x)*log(x)*log(1d0+x)**2)/2d0 - (2*dlog(2d0)
     &*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ log(1d0+x)**4/3d0 + (7*log(x)*zeta3)/8d0 
     &- (3*log(1d0+x)*zeta3)/2d0

      case(48)                  !01-11

         ris = pi**4/72d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &- (basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0))
     &+ 6*cli4(dcmplx(0.5d0,0d0)) + 3*basis3(x) 
     &- basis2(x)/2d0 - (3*basis1(x))/2d0 - 2*basis5(x) 
     &+ 3*basis12(x) - (3*basis15(x))/4d0 
     &- 4*basis4(x) - 4*basis6(x) + 6*basis11(x) 
     &- 6*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- basis3_2(x)*log(1d0-x) + basis3_1(x)*log(1d0-x) 
     &+ 3*basis3_8(x)*log(1d0-x) - basis3_4(x)
     &*log(1d0-x) + basis3_7(x)*log(1d0-x) 
     &+ basis3_5(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 + basis2_1(x)*dlog(2d0)*log(1d0-x)
     &-(dlog(2d0)**3
     &*log(1d0-x))/6d0 - (17*pi**2*log(1d0-x)**2)/48d0 
     &- (basis2_3(x)*log(1d0-x)**2)/2d0 
     &+ (basis2_2(x)*log(1d0-x)**2)/2d0 
     &-(basis2_1(x)*log(1d0-x)**2)/2d0 
     &- (dlog(2d0)**2*log(1d0-x)**2)/4d0 + dlog(2d0)*log(1d0-x)**3 
     &- (47*log(1d0-x)**4)/96d0 + (11*log(1d0-x)**3*log(x))/12d0 
     &+ 2*basis3_1(x)*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/2d0 
     &+ dlog(2d0)**3*log(1d0+x) - (pi**2*log(1d0-x)*log(1d0+x))/24d0 
     &- basis2_1(x)*log(1d0-x)*log(1d0+x) + (dlog(2d0)**2*log(1d0-x)
     &*log(1d0+x))/2d0 - dlog(2d0)*log(1d0-x)**2*log(1d0+x) 
     &+ (log(1d0-x)**3*log(1d0+x))/24d0 - (5*log(1d0-x)**2*log(x)
     &*log(1d0+x))/4d0 + (29*pi**2*log(1d0+x)**2)/48d0 
     &- (3*dlog(2d0)**2*log(1d0+x)**2)/2d0 + (9*log(1d0-x)**2
     &*log(1d0+x)**2)/16d0 - (log(1d0-x)*log(x)*log(1d0+x)**2)/4d0 
     &+ dlog(2d0)*log(1d0+x)**3 + (log(1d0-x)*log(1d0+x)**3)/24d0 
     &+ (7*log(x)*log(1d0+x)**3)/12d0 - (55*log(1d0+x)**4)/96d0 
     &- (29*log(1d0-x)*zeta3)/8d0 + (3*log(1d0+x)*zeta3)/2d0

      case(49)                  !010-1

         ris = basis16(x)

      case(50)                  !0100

         ris = 3*basis1(x) - 2*basis3_1(x)*log(x) 
     &+ (basis2_1(x)*log(x)**2)/2d0

      case(51)                  !0101

         ris = -pi**4/45d0 + basis2_1(x)**2/2d0 + 2*basis3(x) 
     &- 2*basis1(x) - 2*basis5(x) + 2*basis3_1(x)*log(1d0-x) 
     &- (pi**2*log(1d0-x)**2)/6d0 - log(1d0-x)**4/12d0 
     &+ (log(1d0-x)**3*log(x))/3d0 - 2*log(1d0-x)*zeta3

      case(52)                  !011-1

         ris = basis17(x)

      case(53)                  !0110

         ris = -basis2_1(x)**2/2d0 - basis3_3(x)*log(x) 
     &+ (pi**2*log(1d0-x)*log(x))/6d0 - basis2_1(x)*log(1d0-x)*log(x) 
     &- (log(1d0-x)**2*log(x)**2)/2d0 + log(x)*zeta3

      case(54)                  !0111

         ris = pi**4/90d0 - basis3(x) + basis3_3(x)*log(1d0-x) 
     &- (pi**2*log(1d0-x)**2)/12d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (log(1d0-x)**3*log(x))/3d0

      case(55)                  !1-1-1-1

         ris = cli4(dcmplx(0.5d0,0d0)) - basis7(x) 
     &+ basis3_5(x)
     &*log(1d0+x) - (pi**2*log(1d0+x)**2)/12d0 + (basis2_3(x)
     &*log(1d0+x)**2)/2d0 +(dlog(2d0)**2*log(1d0+x)**2)/2d0 - (dlog(2d0)
     &*log(1d0-x)*log(1d0+x)**2)/2d0 - (dlog(2d0)*log(1d0+x)**3)/3d0 
     &+ (log(1d0-x)*log(1d0+x)**3)/3d0

      case(56)                  !1-1-10

         ris = -pi**4/480d0 + basis2_2(x)**2/2d0 
     &+ basis18(x) 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0
     &+ 3*cli4(dcmplx(0.5d0,0d0)) - basis3(x) 
     &- basis2(x) - basis1(x) - basis15(x)/2d0 
     &+ basis9(x)/2d0 - basis10(x)/2d0 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 3*basis7(x) + basis13(x)/4d0 
     &+ basis3_5(x)*log(x) + (pi**2*dlog(2d0)*log(x))/12d0 
     &- (dlog(2d0)**3*log(x))/6d0 + 2*basis3_1(x)*log(1d0+x) 
     &+ 2*basis3_8(x)*log(1d0+x) + 2*basis3_7(x)
     &*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/4d0 + (dlog(2d0)**3
     &*log(1d0+x))/2d0 - (pi**2*log(1d0-x)*log(1d0+x))/4d0 
     &+ dlog(2d0)*log(1d0-x)**2*log(1d0+x) - (log(1d0-x)**3
     &*log(1d0+x))/3d0 - (pi**2*log(x)*log(1d0+x))/6d0 
     &+ basis2_3(x)*log(x)*log(1d0+x) + dlog(2d0)**2*log(x)
     &*log(1d0+x) - dlog(2d0)*log(1d0-x)*log(x)*log(1d0+x) 
     &+ log(1d0-x)**2*log(x)*log(1d0+x) + (5*pi**2*log(1d0+x)**2)
     &/12d0 - (basis2_3(x)*log(1d0+x)**2)/2d0 + (basis2_2(x)
     &*log(1d0+x)**2)/2d0 - (basis2_1(x)*log(1d0+x)**2)/2d0
     &-dlog(2d0)**2
     &*log(1d0+x)**2 - (3*dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 
     &- (dlog(2d0)*log(x)*log(1d0+x)**2)/2d0 - log(1d0-x)*log(x)
     &*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0+x)**3)/2d0 + (log(1d0-x)
     &*log(1d0+x)**3)/2d0 + log(x)*log(1d0+x)**3 - (5*log(1d0+x)**4)
     &/8d0 - (log(1d0-x)*zeta3)/8d0 - (7*log(x)*zeta3)/8d0 
     &+ (5*log(1d0+x)*zeta3)/4d0

      case(57)                  !1-1-11

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*dlog(2d0)**2)/24d0 
     &- (basis2_3(x)*dlog(2d0)**2)/2d0 - dlog(2d0)**4/8d0 
     &- basis3_5(x)*log(1d0-x) - (pi**2*dlog(2d0)*log(1d0-x))/6d0 
     &+ basis2_3(x)*dlog(2d0)*log(1d0-x) + (2*dlog(2d0)**3
     &*log(1d0-x))/3d0 - (dlog(2d0)**2*log(1d0-x)**2)/2d0 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/6d0 - basis2_3(x)
     &*log(1d0-x)*log(1d0+x) - dlog(2d0)**2*log(1d0-x)*log(1d0+x) 
     &+ dlog(2d0)*log(1d0-x)**2*log(1d0+x) + (dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 - (log(1d0-x)**2*log(1d0+x)**2)/2d0 
     &+ (7*log(1d0-x)*zeta3)/8d0

      case(58)                  !1-10-1

         ris = (11*pi**4)/720d0 - basis2_2(x)**2/2d0 
     &- basis18(x) 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*dlog(2d0)**2)/2d0 + basis3_4(x)*log(1d0-x) - basis2_2(x)
     &*dlog(2d0)*log(1d0-x) - 2*basis3_1(x)*log(1d0+x) 
     &- 2*basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) + (pi**2*dlog(2d0)*log(1d0+x))/4d0 - (dlog(2d0)**3
     &*log(1d0+x))/2d0 + (5*pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &+ basis2_2(x)*log(1d0-x)*log(1d0+x) - dlog(2d0)*log(1d0-x)**2
     &*log(1d0+x) + (log(1d0-x)**3*log(1d0+x))/3d0 - log(1d0-x)**2
     &*log(x)*log(1d0+x) - (5*pi**2*log(1d0+x)**2)/12d0 
     &+ (basis2_3(x)*log(1d0+x)**2)/2d0 - (basis2_2(x)
     &*log(1d0+x)**2)/2d0 + (basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ dlog(2d0)**2*log(1d0+x)**2 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + 2*log(1d0-x)*log(x)*log(1d0+x)**2 
     &- (3*dlog(2d0)*log(1d0+x)**3)/2d0 - (2*log(1d0-x)
     &*log(1d0+x)**3)/3d0 - log(x)*log(1d0+x)**3 
     &+ (5*log(1d0+x)**4)/8d0 - (7*log(1d0-x)*zeta3)/8d0 
     &- (5*log(1d0+x)*zeta3)/4d0)
     &+ 2*basis3(x) + (3*basis2(x))/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 - basis4(x) 
     &- basis9(x) + basis10(x) 
     &+ 2*basis6(x) - basis11(x) - basis13(x)/2d0 
     &- basis3_6(x)*log(1d0+x) + basis3_3(x)*log(1d0+x) 
     &- basis3_2(x)*log(1d0+x) - basis3_1(x)*log(1d0+x) 
     &- basis3_8(x)*log(1d0+x) - 2*basis3_7(x)
     &*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/12d0 + (dlog(2d0)**3
     &*log(1d0+x))/6d0 + (pi**2*log(1d0-x)*log(1d0+x))/12d0 
     &- (dlog(2d0)**2*log(1d0-x)*log(1d0+x))/2d0 - (dlog(2d0)
     &*log(1d0-x)**2*log(1d0+x))/2d0 + (log(1d0-x)**3*log(1d0+x))
     &/6d0 - (log(1d0-x)**2*log(x)*log(1d0+x))/2d0 - (pi**2
     &*log(1d0+x)**2)/8d0 + (basis2_3(x)*log(1d0+x)**2)/2d0 
     &- (basis2_2(x)*log(1d0+x)**2)/2d0 
     &+(basis2_1(x)*log(1d0+x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0+x)**2)/4d0 + (3*dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2)/2d0 + log(1d0-x)*log(x)*log(1d0+x)**2 
     &- dlog(2d0)*log(1d0+x)**3 - (log(1d0-x)*log(1d0+x)**3)/2d0 
     &- (5*log(x)*log(1d0+x)**3)/6d0 + (11*log(1d0+x)**4)/24d0 
     &+ (log(1d0-x)*zeta3)/4d0 - (3*log(1d0+x)*zeta3)/8d0

      case(59)                  !1-100

         ris = pi**4/72d0 + 2*cli4(dcmplx(0.5d0,0d0)) + basis3(x) 
     &+ (3*basis2(x))/2d0 - basis1(x)/2d0 - basis5(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (pi**2*log(1d0-x)**2)/16d0 
     &- log(1d0-x)**4/32d0 - basis3_6(x)*log(x) + basis3_3(x)
     &*log(x) - basis3_2(x)*log(x) + basis3_1(x)*log(x) 
     &+ basis3_8(x)*log(x) - (pi**2*dlog(2d0)*log(x))/12d0 
     &+ (dlog(2d0)**3*log(x))/6d0 - (pi**2*log(1d0-x)*log(x))/12d0 
     &- (dlog(2d0)**2*log(1d0-x)*log(x))/2d0 + (dlog(2d0)*log(1d0-x)**2
     &*log(x))/2d0 - (log(1d0-x)**3*log(x))/12d0 - (pi**2*log(x)**2)
     &/24d0 + (basis2_3(x)*log(x)**2)/2d0 + (dlog(2d0)**2
     &*log(x)**2)/4d0 - (dlog(2d0)*log(1d0-x)*log(x)**2)/2d0 
     &+ (log(1d0-x)**2*log(x)**2)/2d0 - (pi**2*dlog(2d0)*log(1d0+x))
     &/6d0 + (dlog(2d0)**3*log(1d0+x))/3d0 + (pi**2*log(1d0-x)
     &*log(1d0+x))/24d0 + (log(1d0-x)**3*log(1d0+x))/24d0 
     &- (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 + (13*pi**2
     &*log(1d0+x)**2)/48d0 - (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &+ (log(1d0-x)**2*log(1d0+x)**2)/16d0 - (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 + (dlog(2d0)*log(1d0+x)**3)/3d0 + (log(1d0-x)
     &*log(1d0+x)**3)/24d0 + (log(x)*log(1d0+x)**3)/4d0 
     &- (23*log(1d0+x)**4)/96d0 - (3*log(1d0-x)*zeta3)/4d0 
     &- (log(x)*zeta3)/8d0

      case(60)                  !1-101

         ris = -pi**4/72d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &-(basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0))
     &+ 4*cli4(dcmplx(0.5d0,0d0)) 
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 + 2*basis11(x) 
     &- 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*log(1d0-x) 
     &- basis3_3(x)*log(1d0-x) + basis3_2(x)*log(1d0-x) 
     &- basis3_1(x)
     &*log(1d0-x) + basis3_8(x)*log(1d0-x) - (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 + (dlog(2d0)**3*log(1d0-x))/6d0 + (5*pi**2
     &*log(1d0-x)**2)/48d0 - (basis2_3(x)*log(1d0-x)**2)/2d0 
     &+ (basis2_2(x)*log(1d0-x)**2)/2d0 
     &-(basis2_1(x)*log(1d0-x)**2)/2d0 
     &- (dlog(2d0)**2*log(1d0-x)**2)/4d0 + (dlog(2d0)*log(1d0-x)**3)/3d0 
     &- (3*log(1d0-x)**4)/32d0 - (log(1d0-x)**3*log(x))/4d0 
     &+ 2*basis3_1(x)*log(1d0+x) - (pi**2*dlog(2d0)*log(1d0+x))/6d0 
     &+ (dlog(2d0)**3*log(1d0+x))/3d0 - (pi**2*log(1d0-x)*log(1d0+x))
     &/24d0 - (log(1d0-x)**3*log(1d0+x))/24d0 + (log(1d0-x)**2
     &*log(x)*log(1d0+x))/4d0 + (pi**2*log(1d0+x)**2)/16d0 
     &- (dlog(2d0)**2*log(1d0+x)**2)/2d0 - (log(1d0-x)**2
     &*log(1d0+x)**2)/16d0 + (log(1d0-x)*log(x)*log(1d0+x)**2)/4d0 
     &+ (dlog(2d0)*log(1d0+x)**3)/3d0 - (log(1d0-x)*log(1d0+x)**3)/24d0 
     &+ (log(x)*log(1d0+x)**3)/12d0 - (3*log(1d0+x)**4)/32d0 
     &+ (5*log(1d0-x)*zeta3)/8d0 + (3*log(1d0+x)*zeta3)/2d0

      case(61)                  !1-11-1

         ris = (-23*pi**4)/1440d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 - 2*basis8(x) 
     &- 2*basis10(x) + 2*basis7(x) 
     &- (pi**2*dlog(2d0)**2)/24d0 + (basis2_3(x)*dlog(2d0)**2)/2d0 
     &+ dlog(2d0)**4/8d0 - (pi**2*dlog(2d0)*log(1d0-x))/12d0 
     &- basis2_3(x)*dlog(2d0)*log(1d0-x) - (dlog(2d0)**3
     &*log(1d0-x))/6d0 + 2*basis3_6(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &+ dlog(2d0)**2*log(1d0-x)*log(1d0+x) - (pi**2*log(1d0+x)**2)/6d0 
     &+ (dlog(2d0)**2*log(1d0+x)**2)/2d0 - dlog(2d0)*log(1d0-x)
     &*log(1d0+x)**2 + (log(1d0-x)*log(1d0+x)**3)/3d0 
     &- log(1d0+x)**4/12d0 + (log(1d0-x)*zeta3)/4d0 
     &- 2*log(1d0+x)*zeta3

      case(62)                  !1-110

         ris = -pi**4/72d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0) - 6*cli4(dcmplx(0.5d0,0d0)) 
     &- 3*basis3(x) + basis2(x)/2d0 + (3*basis1(x))/2d0 
     &+ 2*basis5(x) - 3*basis12(x) 
     &+ (3*basis15(x))/4d0 + 4*basis4(x) 
     &+ 4*basis6(x) - 6*basis11(x) 
     &+ 6*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 - 2*basis3_8(x)
     &*log(1d0-x) + (3*pi**2*log(1d0-x)**2)/16d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/4d0 - (dlog(2d0)*log(1d0-x)**3)/2d0 
     &+ (31*log(1d0-x)**4)/96d0 + 2*basis3_6(x)*log(x) 
     &+ (pi**2*dlog(2d0)*log(x))/6d0 - (dlog(2d0)**3*log(x))/3d0 
     &- (pi**2*log(1d0-x)*log(x))/12d0 - basis2_3(x)
     &*log(1d0-x)*log(x) + (dlog(2d0)**2*log(1d0-x)*log(x))/2d0 
     &- (5*log(1d0-x)**3*log(x))/12d0 - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/2d0 - dlog(2d0)**3*log(1d0+x) 
     &- (pi**2*log(1d0-x)*log(1d0+x))/24d0 - (log(1d0-x)**3
     &*log(1d0+x))/24d0 + (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &- (29*pi**2*log(1d0+x)**2)/48d0 + (3*dlog(2d0)**2
     &*log(1d0+x)**2)/2d0 - (log(1d0-x)**2*log(1d0+x)**2)/16d0 
     &+ (log(1d0-x)*log(x)*log(1d0+x)**2)/4d0 - dlog(2d0)
     &*log(1d0+x)**3 - (log(1d0-x)*log(1d0+x)**3)/24d0 
     &- (7*log(x)*log(1d0+x)**3)/12d0 + (55*log(1d0+x)**4)/96d0 
     &+ (11*log(1d0-x)*zeta3)/4d0 - (7*log(x)*zeta3)/4d0 
     &- (3*log(1d0+x)*zeta3)/2d0

      case(63)                  !1-111

         ris = -3*cli4(dcmplx(0.5d0,0d0)) + 3*basis8(x) 
     &- 2*basis3_6(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - (dlog(2d0)**3*log(1d0-x))/6d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 
     &- (7*log(1d0-x)*zeta3)/8d0

      case(64)                  !10-1-1

         ris = -pi**4/480d0 - basis3(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 + basis13(x)/4d0 
     &+ basis3_6(x)*log(1d0+x) - basis3_3(x)*log(1d0+x) 
     &- basis3_4(x)*log(1d0+x) + basis3_7(x)
     &*log(1d0+x) + basis3_5(x)*log(1d0+x) + (pi**2*dlog(2d0)
     &*log(1d0+x))/6d0 - (dlog(2d0)**3*log(1d0+x))/3d0 + (dlog(2d0)**2
     &*log(1d0-x)*log(1d0+x))/2d0 - (pi**2*log(1d0+x)**2)/6d0 
     &- (basis2_1(x)*log(1d0+x)**2)/2d0
     &+(dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- dlog(2d0)*log(1d0-x)*log(1d0+x)**2 - (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/2d0 + (log(1d0-x)*log(1d0+x)**3)/2d0 
     &- (log(1d0-x)*zeta3)/8d0 - (log(1d0+x)*zeta3)/8d0

      case(65)                  !10-10

         ris = -pi**4/180d0 + basis16(x) - 2*basis3(x) 
     &- 2*basis2(x) - 2*basis1(x) + 2*basis5(x) 
     &+ 2*basis4(x) + 2*basis6(x) + basis13(x)/2d0 
     &- basis14(x)/2d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ log(1d0-x)**4/16d0 + basis3_6(x)*log(x) - basis3_3(x)
     &*log(x) - basis3_4(x)*log(x) + basis3_7(x)
     &*log(x) + basis3_5(x)*log(x) + (pi**2*dlog(2d0)*log(x))
     &/6d0 - (dlog(2d0)**3*log(x))/3d0 - (pi**2*log(1d0-x)*log(x))/12d0 
     &+ (dlog(2d0)**2*log(1d0-x)*log(x))/2d0 - (log(1d0-x)**3
     &*log(x))/6d0 + 2*basis3_1(x)*log(1d0+x) - (pi**2*log(1d0-x)
     &*log(1d0+x))/12d0 - (log(1d0-x)**3*log(1d0+x))/12d0 
     &- (pi**2*log(x)*log(1d0+x))/12d0 - basis2_1(x)*log(x)*log(1d0+x) 
     &+ (dlog(2d0)**2*log(x)*log(1d0+x))/2d0 - dlog(2d0)*log(1d0-x)
     &*log(x)*log(1d0+x) + (log(1d0-x)**2*log(x)*log(1d0+x))/2d0 
     &- log(1d0-x)*log(x)**2*log(1d0+x) - (5*pi**2*log(1d0+x)**2)
     &/24d0 - (log(1d0-x)**2*log(1d0+x)**2)/8d0 + log(1d0-x)*log(x)
     &*log(1d0+x)**2 - (log(1d0-x)*log(1d0+x)**3)/12d0 - (log(x)
     &*log(1d0+x)**3)/6d0 + (7*log(1d0+x)**4)/48d0 + (3*log(1d0-x)
     &*zeta3)/2d0 - (3*log(x)*zeta3)/4d0 + (3*log(1d0+x)*zeta3)/2d0

      case(66)                  !10-11

         ris = pi**4/72d0 - (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &+ basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0)
     &- 4*cli4(dcmplx(0.5d0,0d0)) + 2*basis8(x)
     & - basis3(x) + (3*basis2(x))/2d0 + basis1(x)/2d0 
     &- 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 - 2*basis11(x) 
     &+ 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- basis3_6(x)*log(1d0-x) + basis3_3(x)*log(1d0-x) 
     &- 2*basis3_8(x)*log(1d0-x) + basis3_4(x)
     &*log(1d0-x) - basis3_7(x)*log(1d0-x) 
     &- basis3_5(x)*log(1d0-x) - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &+ (pi**2*log(1d0-x)**2)/16d0 + (basis2_3(x)
     &*log(1d0-x)**2)/2d0 - (basis2_2(x)*log(1d0-x)**2)/2d0 
     &+ (basis2_1(x)*log(1d0-x)**2)/2d0
     &+(dlog(2d0)**2*log(1d0-x)**2)/4d0 
     &- (5*dlog(2d0)*log(1d0-x)**3)/6d0 + (25*log(1d0-x)**4)/96d0 
     &- (log(1d0-x)**3*log(x))/4d0 - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/6d0 - (dlog(2d0)**3*log(1d0+x))/3d0 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/8d0 + basis2_1(x)*log(1d0-x)
     &*log(1d0+x) - (dlog(2d0)**2*log(1d0-x)*log(1d0+x))/2d0 
     &+ dlog(2d0)*log(1d0-x)**2*log(1d0+x) + (log(1d0-x)**3
     &*log(1d0+x))/24d0 + (3*log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &- (pi**2*log(1d0+x)**2)/16d0 + (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- (7*log(1d0-x)**2*log(1d0+x)**2)/16d0 - (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 - (dlog(2d0)*log(1d0+x)**3)/3d0 
     &+ (log(1d0-x)*log(1d0+x)**3)/24d0 - (log(x)*log(1d0+x)**3)
     &/12d0 + (3*log(1d0+x)**4)/32d0 + (log(1d0-x)*zeta3)/4d0 
     &- (3*log(1d0+x)*zeta3)/2d0

      case(67)                  !100-1

         ris = pi**4/360d0 - basis16(x) + basis3(x) + basis2(x) 
     &+ basis1(x) - basis5(x) - basis4(x) 
     &- basis6(x) - basis13(x)/4d0 + basis14(x)
     &/4d0 - (pi**2*log(1d0-x)**2)/16d0 - log(1d0-x)**4/32d0 
     &+ (log(1d0-x)**3*log(x))/12d0 - basis3_1(x)*log(1d0+x) 
     &+ (pi**2*log(1d0-x)*log(1d0+x))/24d0 + (log(1d0-x)**3
     &*log(1d0+x))/24d0 - (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &+ (5*pi**2*log(1d0+x)**2)/48d0 + (log(1d0-x)**2*log(1d0+x)**2)
     &/16d0 - (log(1d0-x)*log(x)*log(1d0+x)**2)/4d0 + (log(1d0-x)
     &*log(1d0+x)**3)/24d0 + (log(x)*log(1d0+x)**3)/12d0 
     &- (7*log(1d0+x)**4)/96d0 - (3*log(1d0-x)*zeta3)/4d0 
     &- (3*log(1d0+x)*zeta3)/4d0

      case(68)                  !1000

         ris = -basis1(x)+basis3_1(x)*log(x)
     &-(basis2_1(x)*log(x)**2)/2d0 
     &- (log(1d0-x)*log(x)**3)/6d0

      case(69)                  !1001

         ris = -basis2_1(x)**2/2d0 - basis3_1(x)*log(1d0-x)

      case(70)                  !101-1

         ris = -pi**4/144d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &+ 4*cli4(dcmplx(0.5d0,0d0)) 
     &- (basis17(x)-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4(dcmplx(0.5d0,0d0)) 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*dlog(2d0)**2)/2d0 
     &- 2*basis3_8(x)*log(1d0-x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 + (pi**2*log(1d0-x)**2)/8d0 
     &+ (basis2_3(x)*log(1d0-x)**2)/2d0 - (basis2_2(x)
     &*log(1d0-x)**2)/2d0 + (basis2_1(x)*log(1d0-x)**2)/2d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/2d0-(2*dlog(2d0)*log(1d0-x)**3)/3d0 
     &+ (7*log(1d0-x)**4)/24d0 - (log(1d0-x)**3*log(x))/3d0 
     &- basis3_3(x)*log(1d0+x) - 2*basis3_1(x)*log(1d0+x) 
     &+ (pi**2*dlog(2d0)*log(1d0+x))/3d0-(2*dlog(2d0)**3*log(1d0+x))/3d0 
     &- (3*pi**2*log(1d0+x)**2)/8d0 + dlog(2d0)**2*log(1d0+x)**2 
     &- (2*dlog(2d0)*log(1d0+x)**3)/3d0 - (log(x)*log(1d0+x)**3)/3d0 
     &+ (3*log(1d0+x)**4)/8d0 + (7*log(1d0-x)*zeta3)/4d0 
     &- (5*log(1d0+x)*zeta3)/8d0))
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &+ basis9(x) - basis10(x) 
     &+ 2*basis11(x) - 2*basis7(x) 
     &- basis13(x)/4d0 - basis14(x)/4d0 
     &- (basis2_1(x)*dlog(2d0)**2)/2d0 + 2*basis3_8(x)*log(1d0-x) 
     &-(pi**2*dlog(2d0)*log(1d0-x))/6d0
     &+basis2_1(x)*dlog(2d0)*log(1d0-x) 
     &+ (dlog(2d0)**3*log(1d0-x))/3d0 - (pi**2*log(1d0-x)**2)/16d0 
     &- (basis2_3(x)*log(1d0-x)**2)/2d0 + (basis2_2(x)
     &*log(1d0-x)**2)/2d0 - (basis2_1(x)*log(1d0-x)**2)/2d0 
     &- (3*dlog(2d0)**2*log(1d0-x)**2)/4d0 + (5*dlog(2d0)*log(1d0-x)**3)
     &/6d0 - (25*log(1d0-x)**4)/96d0 + (log(1d0-x)**3*log(x))/4d0 
     &+ 2*basis3_3(x)*log(1d0+x) + 2*basis3_1(x)*log(1d0+x) - (pi**2
     &*dlog(2d0)*log(1d0+x))/6d0 + (dlog(2d0)**3*log(1d0+x))/3d0 
     &- (pi**2*log(1d0-x)*log(1d0+x))/24d0 - (log(1d0-x)**3
     &*log(1d0+x))/24d0 + (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &+ (7*pi**2*log(1d0+x)**2)/48d0 - (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &- (log(1d0-x)**2*log(1d0+x)**2)/16d0 + (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 + (dlog(2d0)*log(1d0+x)**3)/3d0 - (log(1d0-x)
     &*log(1d0+x)**3)/24d0 + (log(x)*log(1d0+x)**3)/12d0 
     &- (17*log(1d0+x)**4)/96d0 - (3*log(1d0-x)*zeta3)/4d0 
     &- (log(1d0+x)*zeta3)/4d0

      case(71)                  !1010

         ris = pi**4/45d0 + basis2_1(x)**2/2d0 - 2*basis3(x) 
     &+ 2*basis1(x) + 2*basis5(x) + (pi**2*log(1d0-x)**2)/6d0 
     &+ log(1d0-x)**4/12d0 + 2*basis3_3(x)*log(x) - (pi**2
     &*log(1d0-x)*log(x))/3d0 + basis2_1(x)*log(1d0-x)*log(x) 
     &- (log(1d0-x)**3*log(x))/3d0 + log(1d0-x)**2*log(x)**2 
     &+ 2*log(1d0-x)*zeta3 - 2*log(x)*zeta3

      case(72)                  !1011

         ris = -pi**4/30d0 + 3*basis3(x) - 2*basis3_3(x)*log(1d0-x) 
     &+ (pi**2*log(1d0-x)**2)/12d0 - (basis2_1(x)*log(1d0-x)**2)/2d0 
     &- (log(1d0-x)**3*log(x))/2d0 - log(1d0-x)*zeta3

      case(73)                  !11-1-1

         ris = (7*pi**4)/720d0 + basis8(x) 
     &+ basis10(x) - basis7(x) + (pi**2*dlog(2d0)
     &*log(1d0-x))/12d0 - (dlog(2d0)**3*log(1d0-x))/6d0 + (dlog(2d0)**2
     &*log(1d0-x)**2)/4d0 - basis3_6(x)*log(1d0+x) 
     &- (pi**2*dlog(2d0)*log(1d0+x))/6d0 + (dlog(2d0)**3*log(1d0+x))/3d0 
     &- (dlog(2d0)**2*log(1d0-x)*log(1d0+x))/2d0 + (pi**2
     &*log(1d0+x)**2)/12d0 - (dlog(2d0)**2*log(1d0+x)**2)/4d0 
     &+ (dlog(2d0)*log(1d0-x)*log(1d0+x)**2)/2d0 - (log(1d0-x)
     &*log(1d0+x)**3)/6d0 + log(1d0+x)**4/24d0 - (log(1d0-x)*zeta3)
     &/8d0 + log(1d0+x)*zeta3

      case(74)                  !11-10

         ris = pi**4/72d0 + cli4(dcmplx(0.5d0,0d0)) + basis8(x) 
     &+ basis3(x) + basis2(x)/2d0 - basis1(x)/2d0 
     &- 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (pi**2*dlog(2d0)*log(1d0-x))/12d0 
     &- (dlog(2d0)**3*log(1d0-x))/6d0 - (5*pi**2*log(1d0-x)**2)/48d0 
     &+ (dlog(2d0)**2*log(1d0-x)**2)/4d0 - (dlog(2d0)*log(1d0-x)**3)/6d0 
     &- log(1d0-x)**4/32d0 - basis3_6(x)*log(x) - (pi**2
     &*dlog(2d0)*log(x))/12d0 + (dlog(2d0)**3*log(x))/6d0 + (pi**2
     &*log(1d0-x)*log(x))/12d0 - (dlog(2d0)**2*log(1d0-x)*log(x))/2d0 
     &+ (dlog(2d0)*log(1d0-x)**2*log(x))/2d0 + (log(1d0-x)**3*log(x))
     &/12d0 - (pi**2*dlog(2d0)*log(1d0+x))/6d0+(dlog(2d0)**3*log(1d0+x))
     &/3d0 + (pi**2*log(1d0-x)*log(1d0+x))/24d0 + (log(1d0-x)**3
     &*log(1d0+x))/24d0 - (log(1d0-x)**2*log(x)*log(1d0+x))/4d0 
     &+ (13*pi**2*log(1d0+x)**2)/48d0 - (dlog(2d0)**2*log(1d0+x)**2)/2d0 
     &+ (log(1d0-x)**2*log(1d0+x)**2)/16d0 - (log(1d0-x)*log(x)
     &*log(1d0+x)**2)/4d0 + (dlog(2d0)*log(1d0+x)**3)/3d0 + (log(1d0-x)
     &*log(1d0+x)**3)/24d0 + (log(x)*log(1d0+x)**3)/4d0 
     &- (23*log(1d0+x)**4)/96d0 - (13*log(1d0-x)*zeta3)/8d0 
     &+ (7*log(x)*zeta3)/8d0

      case(75)                  !11-11
         
         ris = 3*cli4(dcmplx(0.5d0,0d0)) - 3*basis8(x) 
     &+ basis3_6(x)*log(1d0-x) - (pi**2*dlog(2d0)*log(1d0-x))/6d0 
     &+ (dlog(2d0)**3*log(1d0-x))/3d0 + (pi**2*log(1d0-x)**2)/24d0 
     &- (dlog(2d0)**2*log(1d0-x)**2)/4d0 + (7*log(1d0-x)*zeta3)/4d0

      case(76)                  !110-1

         ris = -pi**4/288d0 + basis4(x) 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ basis13(x)/4d0 + (pi**2*log(1d0-x)**2)/24d0 
     &- basis3_3(x)*log(1d0+x) - (pi**2*log(1d0+x)**2)/24d0 
     &+ log(1d0+x)**4/24d0 + (5*log(1d0-x)*zeta3)/8d0 
     &+ (7*log(1d0+x)*zeta3)/8d0

      case(77)                  !1100

         ris = -pi**4/90d0 + basis3(x) - basis1(x) - basis5(x) 
     &- (pi**2*log(1d0-x)**2)/12d0 - log(1d0-x)**4/24d0 - basis3_3(x)
     &*log(x) + (pi**2*log(1d0-x)*log(x))/6d0 + (log(1d0-x)**3
     &*log(x))/6d0 - (log(1d0-x)**2*log(x)**2)/4d0 
     &- log(1d0-x)*zeta3 + log(x)*zeta3

      case(78)                  !1101

         ris = pi**4/30d0 - 3*basis3(x) + basis3_3(x)*log(1d0-x) 
     &+ (pi**2*log(1d0-x)**2)/12d0 + 2*log(1d0-x)*zeta3

      case(79)                  !111-1

         ris = -cli4(dcmplx(0.5d0,0d0)) + basis8(x) 
     &+ (pi**2*dlog(2d0)*log(1d0-x))/12d0 -(dlog(2d0)**3*log(1d0-x))/6d0 
     &- (pi**2*log(1d0-x)**2)/24d0 + (dlog(2d0)**2*log(1d0-x)**2)/4d0 
     &- (dlog(2d0)*log(1d0-x)**3)/6d0 - (7*log(1d0-x)*zeta3)/8d0

      case(80)                  !1110

         ris = -pi**4/90d0 + basis3(x) - (pi**2*log(1d0-x)**2)/12d0 
     &- log(1d0-x)*zeta3

      case(81)                  !1111

         ris = log(1d0-x)**4/24d0
         
      end select
      endif
      HPL4=ris
      !print*,"HPL4(",n1,",",n2,",",n3,",",n4,",",x,")=",HPL4
      return
      end

c     #####################################################
c     basis1(x) = cli4(x) 
c     basis2(x) = cli4(-x)
c     basis3(x) = cli4(1-x)
c     basis4(x) = cli4(1/(1+x)) 
c     basis5(x) = cli4(x/(x-1))
c     basis6(x) = cli4(x/(x+1)) 
c     basis7(x) = cli4((1+x)/2) 
c     basis8(x) = cli4((1-x)/2)
c     basis9(x) = cli4((1-x)/(1+x))
c     basis10(x) = cli4((x-1)/(x+1))
c     basis11(x) = cli4(2x/(1+x))
c     basis12(x) = cli4(2x/(x-1)) 
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
c     basis14(x) = cli4(x^2/(x^2-1)) 
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
c     basis16(x) = ch2m2(x) 
c     basis17(x) = ch21m1(x) 
c     basis18(x) = ch21m1(-x) 
c     #####################################################

c ---------------------------------------------------------
c     basis1(x) = cli4(x) 
      double complex function basis1(x)
      implicit none
      double complex x,cli4
      basis1=cli4(x)
      return
      end
c ---------------------------------------------------------
c     basis2(x) = cli4(-x)
      double complex function basis2(x)
      implicit none
      double complex x,cli4
      basis2=cli4(-x)
      return
      end
c ---------------------------------------------------------
c     basis3(x) = cli4(1-x)
      double complex function basis3(x)
      implicit none
      double complex x,cli4
      basis3 = cli4(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis4(x) = cli4(1/(1+x)) 
      double complex function basis4(x)
      implicit none
      double complex x,cli4
      basis4 = cli4(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis5(x) = cli4(x/(x-1))
      double complex function basis5(x)
      implicit none
      double complex x,cli4
      basis5 = cli4(x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis6(x) = cli4(x/(x+1))
      double complex function basis6(x)
      implicit none
      double complex x,cli4
      basis6 = cli4(x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis7(x) = cli4((1+x)/2) 
      double complex function basis7(x)
      implicit none
      double complex x,cli4
      basis7 = cli4((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis8(x) = cli4((1-x)/2)
      double complex function basis8(x)
      implicit none
      double complex x,cli4
      basis8 = cli4((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis9(x) = cli4((1-x)/(1+x))
      double complex function basis9(x)
      implicit none
      double complex x,cli4
      basis9 = cli4((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis10(x) = cli4((x-1)/(x+1))
      double complex function basis10(x)
      implicit none
      double complex x,cli4
      basis10 = cli4((x-1d0)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis11(x) = cli4(2x/(1+x))
      double complex function basis11(x)
      implicit none
      double complex x,cli4
      basis11 = cli4(2d0*x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis12(x) = cli4(2x/(x-1)) 
      double complex function basis12(x)
      implicit none
      double complex x,cli4
      basis12 = cli4(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
      double complex function basis13(x)
      implicit none
      double complex x,cli4_sbc !,cli4
      basis13=cli4_sbc(x)
      return
      end
c ---------------------------------------------------------
c     basis14(x) = cli4(x^2/(x^2-1)) 
      double complex function basis14(x)
      implicit none
      double complex x,cli4
      basis14 = cli4(x**2/(x**2-1d0))
      return
      end
c ---------------------------------------------------------
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
      double complex function basis15(x)
      implicit none      
      double complex x,cli4_sbc_2
      basis15=cli4_sbc_2(x)
      return
      end
c ---------------------------------------------------------
c     basis16(x) = ch2m2(x) 
      double complex function basis16(x)
      implicit none
      double complex x,ch2m2
      basis16=ch2m2(x)
      return
      end
c ---------------------------------------------------------
c     basis17(x) = ch21m1(x) 
      double complex function basis17(x)
      implicit none
      double complex x,ch21m1
      basis17=ch21m1(x)
      return
      end
c ---------------------------------------------------------
c     basis18(x) = ch21m1(-x)
      double complex function basis18(x)
      implicit none
      double complex x,ch21m1
      basis18=ch21m1(-x)
      return
      end
c ---------------------------------------------------------
C-------------------------------------------------------------------------
C---  mapping of the tetralog into the convergent region
      
      double complex  function cli4(z)
      implicit none
      double complex ris, z, bsli4_outside, bsli4_inside, wcli4, myi
      double precision zabs, pi, zeta2, zeta3, zeta4, border,zim
      integer signim
     
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      
      border = 0.3d0
      zabs = abs(z)
      zim = dimag(z)
      if (dimag(z).ne.0d0) then 
         signim = dimag(z)/dabs(dimag(z)) ! the sign of the imaginary part of z (+- 1)
      else
         signim = 1 ! for real values, we use the "+i epsilon" prescription, therefore signim is +1
      endif
      
      if (z.eq.dcmplx(1d0,0d0)) then 
         ris = dcmplx(zeta4,0d0)
      else 
         if (zabs.le.1d0) then
            if (zabs.le.border) then ! on the mini-disc (remember those?) of |z| <= 0.3, we use the log(1-x) expansion
               ris=bsli4_inside(z)
            else ! on the annulus 0.3 < |z| <= 1, we use the log(x) expansion
               ris=bsli4_outside(z)
            endif
         else                   ! outside the unit circle, the inversion mapping is needed. NOTE: this is "our" mapping, derived from the integral rep. of Li4 and absorbing
c                                 all imaginary parts into log(-z).
            ris=-wcli4(1d0/z) -log(-z)**4/24d0 - 7d0*zeta4/4d0 
     &           - zeta2*log(-z)**2/2d0
c     --- claude
c     ris = -wcli4(1d0/z) - 1d0/24d0*log(1d0/z)**4 
c     &           + zeta2*log(1d0/z)**2 - myi*pi*signim/6d0*log(1d0/z)**3
c     &           *(-z/zabs)
c     &           + 2d0*zeta4
         endif
      endif
      
      cli4=ris
      return
      end
      
c     --- recursion for li4
      
      double complex  function wcli4(z)
      implicit none
      double complex z, cli4
      wcli4 = cli4(z)
      return
      end

c --- the case Li4(1-z^2) needs some special treatment because of its branch cut structure
c --- (that's what 'sbc' stands for: special branch cut)

      double complex function cli4_sbc(z)
      implicit none
      double complex ris, z, cli4, myi,basis14 !,cli4_with_signim
      double complex ll1,ll2,ll3
      double precision pi,zabs,zreal
      integer signim
      
      pi=3.1415926535897932385D0
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
      if (dimag(z).ne.0d0) then 
         signim = dimag(z)/dabs(dimag(z)) ! the sign of the imaginary part of z (+- 1)
      else
         signim = 1
      endif
           
      if (zabs.le.1d0) then !normal li4
         if (zreal.gt.0d0) then
            ris = cli4(1d0 - z**2)
         else if (zreal.eq.0d0 .and. signim.eq.1) then !also normal li4
            ris = cli4(1d0 - z**2)
         else                   ! special branch cut configuration
            ris = cli4(1d0 - z**2)
c     &           - myi*pi*signim/3d0*log(1d0 - z**2)**3
     &           - myi*pi*signim/3d0*(log(1d0 - z)+log(1d0+z))**3 !claude's suggestions for stabilisation
         endif
      else !inversion mapping needed ! this is Claude's inversion mapping, featuring an explicit imaginary part
         ! note that this inversion mapping doesnt map into itself again, i.e. in claude's notation: B_4^(13)(x) = (...) - B_4^(14)(1/x)
         ! --> we dont need a recursion here
         ll1=log(1d0/z)
         ll2=log(1d0 - 1d0/z)
         ll3=log(1d0 + 1d0/z)
         ris = -2d0/3d0*ll1**4 + 4d0/3d0*ll2
     &        *ll1**3 + 4d0/3d0*ll3*ll1**3 
     &        - ll2**2*ll1**2 - ll3**2
     &        *ll1**2 - 2d0*ll2*ll3
     &        *ll1**2 - pi**2/3d0*ll1**2 + 1d0/3d0
     &        *ll2**3*ll1 + 1d0/3d0
     &        *ll3**3*ll1 + ll2
     &        *ll3**2*ll1 + pi**2/3d0*ll1
     &        *ll2 + ll3*ll2**2
     &        *ll1 + pi**2/3d0*ll1*ll3 
     &        - 1d0/24d0*ll2**4 - 1d0/24d0
     &        *ll3**4 - 1d0/6d0*ll2
     &        *ll3**3 - pi**2/12d0*ll2**2 
     &        - pi**2/12d0*ll3**2 
     &        - 1d0/4d0*ll2**2*ll3**2 
     &        - 1d0/6d0*ll2**3*ll3 
     &        - pi**2/6d0*ll2*ll3 
     &        - 7*pi**4/360d0 
     &        -basis14(1d0/z)
      endif
      
      cli4_sbc = ris
      return
      end


c --- the case Li4(4z/(1+z)^2) also needs some special treatment because of its branch cut structure

      double complex function cli4_sbc_2(z)
      implicit none
      double complex ris, z, cli4, myi, wcli4_sbc_2
      double complex arg
      double precision pi,zabs,zreal
      integer signim
      
      pi=3.1415926535897932385D0
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
      if (dimag(z).ne.0d0) then 
         signim = dimag(z)/dabs(dimag(z))
      else
         signim = 1
      endif
      
      
      ris = dcmplx(0d0,0d0)
      if (zabs.eq.1d0) then !see p. 19 of claudes notes.
         arg = dcmplx(dreal(4d0*z/(1d0+z)**2),signim*1d-60)
         ris = cli4(arg)
      else
         if (zabs.lt.1d0) then
            ris = cli4(4d0*z/(1d0+z)**2)
         else                   !inversion mapping needed !
            ris = wcli4_sbc_2(1d0/z) + myi*pi*signim*
     &(4d0*dlog(2d0)**2*log(1d0/z) - 8d0*dlog(2d0)**2*log(1d0+1d0/z) 
     &+ 2d0*dlog(2d0)*log(1d0/z)**2 - 8d0*dlog(2d0)*log(1d0/z)
     &*log(1d0+1d0/z) + 8d0*dlog(2d0)*log(1d0+1d0/z)**2 
     &+ 1d0/3d0*log(1d0/z)**3 - 2d0*log(1d0+1d0/z)*log(1d0/z)**2 
     &+ 4d0*log(1d0+1d0/z)**2*log(1d0/z) - 8d0/3d0*log(1d0+1d0/z)**3 
     &+ 8d0/3d0*dlog(2d0)**3)
         endif
      endif
      
      cli4_sbc_2 = ris
      return
      end

c     --- recursion for cli4_sbc_2
      
      double complex  function wcli4_sbc_2(z)
      implicit none
      double complex z, cli4_sbc_2
      wcli4_sbc_2 =  cli4_sbc_2(z)
      return
      end

C-------------------------------------------------------------------------
C     mapping of H_2-2(z) into convergent region
      
      double complex  function ch2m2(z)
      implicit none
      double complex ris,z,bsh2m2_inside,bsh2m2_outside,cli4,xcli2
      double complex HPL4,wch2m2,myi !,cli4_sbc,cli3
      double precision pi,zeta2,zeta3,zeta4,zabs,zreal,border
      integer signim

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      
      border = 0.3d0
      zabs = abs(z)
      zreal = dreal(z)
      if (dimag(z).ne.0d0) then 
         signim = dimag(z)/dabs(dimag(z)) ! the sign of the imaginary part of z (+- 1)
      else
         signim = 1
      endif
      
      
      if (zabs.le.1d0) then
         if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
            ris = bsh2m2_inside(z)
         else
            if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
               ris = bsh2m2_outside(z)
            else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
               ris = HPL4(0,-1,0,1,-z) 
            endif
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
c     NOTE: inversion formula from HPL with additional sign(im(z)) dependence
         ris = wch2m2(1d0/z) + 37d0*pi**4/720d0 
     &        - HPL4(0,1,0,0,1d0/z) 
     &        - log(1d0/z)**4/24d0 - pi**2/12d0*log(1d0/z)**2 
     &        - pi**2/6d0*xcli2(1d0/z) 
     &        - cli4(-1d0/z) 
     &        + 3d0*zeta3*log(1d0/z)/2d0 
     &        - pi**3*myi*signim*log(1d0/z)/12d0
      endif
      ch2m2=ris
      return
      end
      
      
c     --- recursion for H_2-2(z)
      
      double complex  function wch2m2(z)
      implicit none
      double complex z, ch2m2
           
      wch2m2 = ch2m2(z)
      return
      end

C-------------------------------------------------------------------------
C     mapping of H21-1(z) into convergent region 
      
      double complex  function ch21m1(z)
      implicit none
      double complex ris,z,bsh21m1_inside,bsh21m1_outside_1
      double complex bsh21m1_outside_2,cli4,xcli2,HPL4,wch21m1,myi,ch2m2
      double precision pi,zeta2,zeta3,zeta4,border,zreal,zabs
      integer signim

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      border = 0.3d0
      myi = dcmplx(0d0,1d0)

      if (dimag(z).ne.0d0) then 
         signim = dimag(z)/dabs(dimag(z))
      else
         signim = 1
      endif

      zabs = abs(z)
      zreal = dreal(z)
      
      
      if (zabs.le.1d0) then
         if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1-z) expansion
            ris = bsh21m1_inside(z)
         else
            if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have a log(x) exp.
               ris = bsh21m1_outside_1(z)
            else                ! for Re(z) < 0, there is a different log(x) expansion
               ris = bsh21m1_outside_2(z)
            endif
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
c     Note: inversion formula from HPL with additional sign(im(z)) dependence
         ris = -wch21m1(1d0/z)-pi**4/144d0 -ch2m2(1d0/z) 
     &        - HPL4(0,0,1,-1,1d0/z) + HPL4(0,0,1,0,1d0/z) 
     &        + HPL4(0,1,0,0,1d0/z) 
     &        + HPL4(0,1,1,0,1d0/z) + log(1d0/z)**4/24d0 
     &        + pi**2*dlog(2d0)**2/3d0 - dlog(2d0)**4/12d0 
     &        + 3d0*pi**2*dlog(2d0)*log(1d0/z)/4d0 
     &        + pi**2*log(1d0/z)**2/8d0 
     &        + pi**2*xcli2(1d0/z)/4d0 
     &        - 2*cli4(dcmplx(0.5d0,0d0)) 
     &        + cli4(-1d0/z) 
     &        - 7d0*zeta3*log(1d0/z)/8d0 
     &        + myi*signim*(pi**3*dlog(2d0)/6d0 
     &        + pi**3*log(1d0/z)/12d0 
     &        - 0.5d0*pi*dlog(2d0)**2*log(1d0/z) 
     &        - 0.5d0*pi*dlog(2d0)*log(1d0/z)**2 
     &        - pi*dlog(2d0)*xcli2(1d0/z))
         
      endif
      ch21m1=ris
      return
      end
      

c     --- recursion for H21-1(z)
      
      double complex  function wch21m1(z)
      implicit none
      double complex z, ch21m1
            
      wch21m1 =  ch21m1(z)
      return
      end
      
C---- expansion of tetralogarithm in y = - log(1-z) with Bernoulli numbers  
C------requires  routine fbern4 in coefficients.F for the  coefficients 
C-------of the series  expansion 
      
      double  complex function bsli4_inside(z)
      implicit none
      integer i, Nmin
      double complex elem, ris, z, zb 
      double precision coeffi, fbern4
      Nmin=60
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmin 
         coeffi=fbern4(i)
         elem = zb**(i+1)*coeffi/(i+1)
         ris = ris+elem
      enddo
      
      bsli4_inside=ris 
      return 
      end

C---- expansion of tetralogarithm in y = log(z) with Zeta values  
C------requires  routine zetaval in coefficients.F for the  coefficients 
C-------of the series  expansion 
C-------- used for 0.3 < |z| < 1
      
      double  complex function bsli4_outside(z)
      implicit none
      integer i, Nmin
      double complex elem, ris, z, zb, coeffi
      double precision zetaval4
      Nmin=60
      zb = log(z)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmin 
         if (i.eq.3) then
            coeffi = 1d0 + 0.5d0 + 1d0/3d0 -log(-zb)
            coeffi = coeffi/6d0 !factorial, for zeta values absorbed in zetaval
         else
            coeffi = dcmplx(zetaval4(i),0d0)
         endif
         elem = zb**i*coeffi
         ris = ris+elem
      enddo
      
      bsli4_outside=ris 
      return 
      end

C---- expansion of H2m2(z) = -Li22(-1,z) in y=-log(1+z)
C---- requires the routine bsh2m2_inside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh2m2_inside(z)
      implicit none
      double complex z,y,elem,ris,coeff
      double precision bsh2m2_inside_coeff
      integer n, Nmax

      Nmax = 50
      
      y = dcmplx(1d0,0d0)+z
      y = -log(y)
      ris = dcmplx(0d0,0d0)
      do n=0,Nmax
         coeff = bsh2m2_inside_coeff(n)
         elem = y**(n+1)*coeff
         ris = ris + elem
      enddo
      
      bsh2m2_inside = ris
      return 
      end

C---- expansion of H2m2(z) = -Li22(-1,z) in y = log(z) (and Re(z) >= 0)
C---- requires the routine bsh2m2_outside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh2m2_outside(z)
      implicit none
      double complex z,y,elem,ris,coeff,xcli2,cli4
      double precision bsh2m2_outside_coeff,pi,zeta3
      integer n, Nmax

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      Nmax = 50
   
      y = log(z)
      ris = dcmplx(0d0,0d0)
      do n=0,Nmax
         coeff = bsh2m2_outside_coeff(n)
         elem = y**(n+2)*coeff
         ris = ris + elem
      enddo
      ! additional pieces (not part of the sum)
      ris = ris + 71d0/1440d0*pi**4 + 1d0/6d0*pi**2*dlog(2d0)**2 
     &     - dlog(2d0)**4/6d0 - 4d0*cli4(dcmplx(0.5d0,0d0)) 
     &     - 7d0/2d0*dlog(2d0)*zeta3 
     &     - 5d0/8d0*zeta3*log(z) + pi**2/12d0*xcli2(z) - pi**4/72d0
      
      bsh2m2_outside = ris
      return 
      end


C---- expansion of H_2,1,-1(z) in y = log(1-z)
C---- requires the routine bsh21m1_inside_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_inside(z)
      implicit none
      double complex z,y,elem,ris,coeff
      double precision bsh21m1_inside_coeff
      integer n, Nmax
      
      Nmax = 50
      
      y = dcmplx(1d0,0d0)-z
      y = -log(y)
      ris = dcmplx(0d0,0d0)
      do n=0,Nmax
         coeff = bsh21m1_inside_coeff(n)
         elem = y**(n+1)*coeff
         ris = ris + elem
      enddo

      ris = ris 

      bsh21m1_inside = ris
      return
      end

C---- expansion of H_2,1,-1(z) in y = log(z) for Re(z) >= 0
C---- requires the routine bsh21m1_outside_1_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_outside_1(z)
      implicit none
      double complex z,y,elem,ris,coeff,xcli2,cli4,cli3
      double precision bsh21m1_outside_1_coeff,pi,zeta3
      integer n, Nmax

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      Nmax = 50

      y = log(z)
      ris = dcmplx(0d0,0d0)
      do n=0,Nmax
         coeff = bsh21m1_outside_1_coeff(n)
         elem = y**(n+2)*coeff
         ris = ris + elem
      enddo
      
      ! additional pieces (not part of the sum)
      ris = ris - pi**4/80d0 + pi**2/12d0*dlog(2d0)**2+dlog(2d0)**4/24d0 
     &     + cli4(dcmplx(0.5d0,0d0)) + 7d0/8d0*dlog(2d0)*zeta3 
     &     + y*(7d0*zeta3/8d0 + dlog(2d0)**3/6d0 
     &     - pi**2/12d0*dlog(2d0)) - (0.5d0*dlog(2d0)**2 - pi**2/12d0)
     &     *(pi**2/6d0 - xcli2(z)) + dlog(2d0)*(-cli3(1d0 - z) 
     &     + xcli2(1d0 - z)*log(1d0 - z) + 0.5d0*y*log(1d0-z)**2)
   
      bsh21m1_outside_1 = ris
      return
      end

C---- expansion of H_2,1,-1(z) in y = log(-z) for Re(z) <= 0
C---- requires the routine bsh21m1_outside_2_coeff in li22coeff.F for the coefficients
C---- Nmax is the highest order of the taylor expansion

      double complex function bsh21m1_outside_2(z)
      implicit none
      double complex z,y,elem,ris,coeff1,coeff2,coeff3,xcli2,cli4,m1
      double precision bsh21m1_outside_2_coeff,pi,zeta3
      integer n, Nmax

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      Nmax = 50
      m1 = dcmplx(-1d0,0d0)
      
      y = log(-z)
      ris = dcmplx(0d0,0d0)
      
      do n=0,Nmax
         coeff1 = bsh21m1_outside_2_coeff(n,1)
     &        *(log(-y) - 1d0/(n+1) - 1d0/(n+2))/4d0
         coeff2 = -bsh21m1_outside_2_coeff(n,2)/4d0
         coeff3 = -bsh21m1_outside_2_coeff(n,3)/4d0
         elem = y**(n+2)*(coeff1 + coeff2 + coeff3)
         ris = ris + elem
      enddo
      
!     additional pieces (not part of the sum)
      ris = ris + pi**4/80d0 - pi**2*dlog(2d0)**2/24d0 
     &     - dlog(2d0)**4/12d0 - (pi**2/12d0 
     &     - dlog(2d0)**2/2d0)*(-pi**2/12d0 -dlog(2d0)*y - xcli2(z)) 
     &     - 2d0*cli4(dcmplx(0.5d0,0d0)) 
     &     - y*(-dlog(2d0)**3/6d0 + zeta3/8d0)
      
      bsh21m1_outside_2 = ris 
      return
      end
