C---- Evaluation of rank 3 HPLs with complex argument 

c --- Real part of HPL3     
      double precision function HPL3real(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3real = dreal(HPL3(n1,n2,n3,x))
      return
      end

c --- Imaginary part of HPL3     
      double precision function HPL3im(n1,n2,n3,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2,n3
      double complex x,HPL3
      x=dcmplx(xr,xi)
      HPL3im = dimag(HPL3(n1,n2,n3,x))
      return
      end

C---  maps to basis function (trilogarithm)

      double complex function HPL3(n1, n2, n3, x)
      implicit none
      double precision pi, zeta2, zeta3
      double complex x, ris,myi
      double complex basis3_1,basis3_2,basis3_3,basis3_4
      double complex basis3_5,basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
c      double complex cli2,cli3 ! not needed ATM due to mapping to basis functions
      integer n1,n2,n3,j 

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9

      ris=dcmplx(10d70,10d70)

      if (x.eq.(dcmplx(0d0,0d0))) then 
         if (j.ne.14) then             
            ris=dcmplx(0d0,0d0)
         elseif (j.eq.14) then
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL3: "
            print*, "HPL3(",n1,",",n2,",",n3,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      elseif (x.eq.(dcmplx(1d0,0d0))) then  
         if(j.le.18.or.j.eq.23) then
            select case (j)
            case (1)
               ris =dlog(2d0)**3d0/6d0
            case (2)
               ris =-(pi**2d0*dlog(2d0))/12d0 + zeta3/8d0
            case (3)
               ris =-dlog(2d0)**3d0/6d0 + zeta3/8d0
            case (4)
               ris =(pi**2d0*dlog(2d0))/12d0 - zeta3/4d0
            case (5)
               ris =(3d0*zeta3)/4d0
            case (6)
               ris =(pi**2d0*dlog(2d0))/6d0 - (5d0*zeta3)/8d0
            case (7)
               ris =(pi**2d0*dlog(2d0))/12d0 
     &              - dlog(2d0)**3d0/6d0-zeta3/4d0
            case (8)
               ris =(pi**2d0*dlog(2d0))/12d0 - zeta3
            case (9)
               ris =-(pi**2*dlog(2d0))/12d0+dlog(2d0)**3/6d0
     &              +(7d0*zeta3)/8d0
            case (10)
               ris =zeta3/8d0
            case (11)
               ris =(-3d0*zeta3)/2d0
            case (12)
               ris =-(pi**2d0*dlog(2d0))/4d0 + (13d0*zeta3)/8d0
            case (13)
               ris =(3d0*zeta3)/4d0
            case (14)
               ris =0d0
            case (15)
               ris =zeta3
            case (16)
               ris =(pi**2d0*dlog(2d0))/4d0 - zeta3
            case (17)
               ris =-2d0*zeta3
            case (18)
               ris =zeta3
            case (23)
               ris =zeta3
            end select
         else
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL3: "
            print*, "HPL3(",n1,",",n2,",",n3,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      elseif (x.eq.-1d0) then
         if(j.ge.10) then
            select case (j)
            case (10)
               ris = zeta3
            case (11)
               ris = 2*zeta3 - myi*pi*zeta2
            case(12)
               ris = pi**2*dlog(2d0)/4d0 - zeta3
            case(13)
               ris = -zeta3
            case(14)
               ris = -myi*pi*zeta2
            case(15)
               ris = 3*zeta3/4d0
            case(16)
               ris = -pi**2*dlog(2d0)/4d0 + 13d0*zeta3/8d0
            case(17)
               ris = 3*zeta3/2d0 - myi*pi**3/12d0
            case(18)
               ris = zeta3/8d0
            case(19)
               ris = 0.5d0*zeta2*dlog(2d0) - dlog(2d0)**3/6d0 
     &              - 7d0*zeta3/8d0
            case(20)
               ris = 0.5d0*zeta2*dlog(2d0) + myi*pi*(0.5d0*zeta2 
     &              - 0.5d0*dlog(2d0)**2) - zeta3
            case(21)
               ris = -0.5d0*zeta2*dlog(2d0) + dlog(2d0)**3/6d0 
     &              + zeta3/4d0
            case(22)
               ris = zeta2*dlog(2d0) - 5d0*zeta3/8d0
            case(23)
               ris = 0.5d0*myi*pi*zeta2 + pi**2*dlog(2d0)/2d0 
     &              - 3*zeta3/4d0
            case(24)
               ris = 0.5d0*zeta2*dlog(2d0) - zeta3/4d0
            case(25)
               ris = dlog(2d0)**3/6d0 - zeta3/8d0
            case(26)
               ris = -0.5d0*zeta2*dlog(2d0) + zeta3/8d0 
     &              + 0.5d0*myi*pi*dlog(2d0)**2
            case(27)
               ris = -1d0/6d0*dlog(2d0)**3
            end select
         else
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL3: "
            print*, "HPL3(",n1,",",n2,",",n3,",",x
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      else

c ########## hack to get branch cuts right #################
         if (dimag(x).eq.0d0) then
            x = x + dcmplx(0d0,1d-60)
         endif
c ##########################################################

         select case(j)

c     #####################################################
c     basis3_1(z) = cli3(z) 
c     basis3_2(z) = cli3(-z)
c     basis3_3(z) = cli3(1-z)
c     basis3_4(z) = cli3(1/(1+z)) 
c     basis3_5(z) = cli3((1+z)/2) 
c     basis3_6(z) = cli3((1-z)/2) 
c     basis3_7(z) = cli3((1-z)/(1+z)) 
c     basis3_8(z) = cli3(2z/(z-1))
c     basis2_1(z) = cli2(z)
c     basis2_2(z) = cli2(-z)) 
c     basis2_3(z) = cli2((1-z)/2)
c     #####################################################

      case(1)
         ris = log(1d0+x)**3/6d0
      case(2)
         ris = -(pi**2*log(1d0+x))/6d0 + log(1d0+x)**3/6d0 
     &- basis3_4(x) 
     &        + zeta3
      case(3)
         ris = (pi**2*dlog(2d0))/12d0 - dlog(2d0)**3/6d0 
     &- (pi**2*log(1d0+x))/12d0 + (dlog(2d0)**2*log(1d0+x))/2d0 
     &- (dlog(2d0)*log(1d0+x)**2)/2d0 + basis3_5(x) - (7*zeta3)/8d0
      case(4)
         ris = (pi**2*log(1d0+x))/3d0 + log(x)*log(1d0+x)**2 
     &- log(1d0+x)**3/3d0 + log(1d0+x)*basis2_2(x)+2*basis3_4(x)
     &-2*zeta3
      case(5)
         ris = (log(x)**2*log(1d0+x))/2d0+log(x)*basis2_2(x)
     &-basis3_2(x)
      case(6)
         ris = (pi**2*dlog(2d0))/6d0 - dlog(2d0)**3/3d0 
     &- (pi**2*log(1d0-x))/12d0 + (dlog(2d0)**2*log(1d0-x))/2d0 
     &- (pi**2*log(1d0+x))/12d0 + (dlog(2d0)**2*log(1d0+x))/2d0 
     &- dlog(2d0)*log(1d0-x)*log(1d0+x) - log(1d0-x)*log(x)*log(1d0+x) 
     &+(log(1d0-x)*log(1d0+x)**2)/2d0-log(1d0-x)*basis2_2(x)
     &+basis3_6(x) 
     &- basis3_3(x)-basis3_4(x)+ basis3_7(x) + basis3_5(x)
     &-(3*zeta3)/4d0
      case(7)
         ris = -(pi**2*dlog(2d0))/6d0 + dlog(2d0)**3/3d0 
     &+ (pi**2*log(1d0+x))/4d0 - (3*dlog(2d0)**2*log(1d0+x))/2d0 
     &+ dlog(2d0)*log(1d0-x)*log(1d0+x) + dlog(2d0)*log(1d0+x)**2 
     &- log(1d0-x)*log(1d0+x)**2 - log(1d0+x)*basis2_3(x) 
     &- 2*basis3_5(x) + (7*zeta3)/4d0
      case(8)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &+ (pi**2*log(1d0-x))/6d0 - (dlog(2d0)*log(1d0-x)**2)/2d0 
     &+ log(1d0-x)**3/6d0+(pi**2*log(x))/12d0- (dlog(2d0)**2*log(x))/2d0 
     &+ dlog(2d0)*log(1d0-x)*log(x) - (log(1d0-x)**2*log(x))/2d0 
     &+ (pi**2*log(1d0+x))/12d0 - (dlog(2d0)**2*log(1d0+x))/2d0 
     &+ dlog(2d0)*log(1d0-x)*log(1d0+x) - (log(1d0-x)*log(1d0+x)**2)/2d0 
     &- log(x)*basis2_3(x) + basis3_2(x)-basis3_1(x)-basis3_8(x) 
     &+ basis3_4(x) - basis3_7(x) - basis3_5(x) + (7*zeta3)/8d0
      case(9)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &- (dlog(2d0)*log(1d0-x)**2)/2d0 + (log(1d0-x)**2*log(1d0+x))/2d0 
     &+ log(1d0-x)*basis2_3(x) - basis3_6(x) + (7*zeta3)/8d0
      case(10)
         ris = -(pi**2*log(1d0+x))/6d0 - (log(x)*log(1d0+x)**2)/2d0 
     &+ log(1d0+x)**3/6d0 - log(1d0+x)*basis2_2(x) - basis3_4(x) 
     &+ zeta3
      case(11)
         ris = -(log(x)*basis2_2(x)) + 2*basis3_2(x)
      case(12)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &- (pi**2*log(1d0-x))/12d0 - (dlog(2d0)**2*log(1d0-x))/2d0 
     &+ (dlog(2d0)*log(1d0-x)**2)/2d0 - log(1d0-x)**3/6d0 
     &+ (log(1d0-x)**2*log(x))/2d0 + log(1d0-x)*basis2_2(x) 
     &- basis3_6(x) + basis3_3(x) - basis3_2(x) +basis3_1(x)
     &+basis3_8(x) 
     &- zeta3/8d0
      case(13)
         ris = -basis3_2(x)
      case(14)
         ris = log(x)**3/6d0
      case(15)
         ris = basis3_1(x)
      case(16)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &+ (pi**2*log(1d0-x))/6d0 - (dlog(2d0)*log(1d0-x)**2)/2d0 
     &+ log(1d0-x)**3/6d0 - (log(1d0-x)**2*log(x))/2d0 
     &+ (pi**2*log(1d0+x))/12d0 - (dlog(2d0)**2*log(1d0+x))/2d0 
     &+ dlog(2d0)*log(1d0-x)*log(1d0+x) + log(1d0-x)*log(x)*log(1d0+x) 
     &-(log(1d0-x)*log(1d0+x)**2)/2d0+log(1d0+x)*basis2_1(x)
     &+basis3_2(x) 
     &- basis3_1(x) - basis3_8(x) + basis3_4(x) -basis3_7(x)
     &-basis3_5(x) 
     &+ (7*zeta3)/8d0
      case(17)
         ris = log(x)*basis2_1(x) - 2*basis3_1(x)
      case(18)
         ris = (pi**2*log(1d0-x))/6d0 - (log(1d0-x)**2*log(x))/2d0 
     &- log(1d0-x)*basis2_1(x) - basis3_3(x) + zeta3
      case(19)
         ris = (pi**2*dlog(2d0))/12d0 - dlog(2d0)**3/6d0 
     &- (pi**2*log(1d0+x))/6d0 + dlog(2d0)**2*log(1d0+x) 
     &- dlog(2d0)*log(1d0-x)*log(1d0+x) - (dlog(2d0)*log(1d0+x)**2)/2d0 
     &+ (log(1d0-x)*log(1d0+x)**2)/2d0 + log(1d0+x)*basis2_3(x) 
     &+ basis3_5(x) - (7*zeta3)/8d0
      case(20)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &- (pi**2*log(1d0-x))/12d0 - (dlog(2d0)**2*log(1d0-x))/2d0 
     &+ (dlog(2d0)*log(1d0-x)**2)/2d0 - log(1d0-x)**3/6d0 
     &- (pi**2*log(x))/12d0 + (dlog(2d0)**2*log(x))/2d0 
     &- dlog(2d0)*log(1d0-x)*log(x) + (log(1d0-x)**2*log(x))/2d0 
     &+ log(x)*basis2_3(x) - basis3_6(x) + basis3_3(x) 
     &- basis3_2(x) + basis3_1(x) + basis3_8(x) - zeta3/8d0
      case(21)
         ris = (pi**2*dlog(2d0))/6d0 - dlog(2d0)**3/3d0 
     &- (pi**2*log(1d0-x))/12d0 + (dlog(2d0)**2*log(1d0-x))/2d0 
     &- log(1d0-x)*basis2_3(x) + 2*basis3_6(x) - (7*zeta3)/4d0
      case(22)
         ris = (pi**2*dlog(2d0))/6d0 - dlog(2d0)**3/3d0 
     &- (pi**2*log(1d0-x))/12d0 + (dlog(2d0)**2*log(1d0-x))/2d0 
     &- (pi**2*log(1d0+x))/12d0 + (dlog(2d0)**2*log(1d0+x))/2d0 
     &- dlog(2d0)*log(1d0-x)*log(1d0+x) - log(1d0-x)*log(x)*log(1d0+x) 
     &+(log(1d0-x)*log(1d0+x)**2)/2d0-log(1d0+x)*basis2_1(x)
     &+basis3_6(x) 
     &- basis3_3(x) - basis3_4(x)+basis3_7(x)+basis3_5(x)
     &- (3*zeta3)/4d0
      case(23)
         ris =-(log(1d0-x)*log(x)**2)/2d0-log(x)*basis2_1(x)
     &+basis3_1(x)
      case(24)
         ris = -(pi**2*log(1d0-x))/3d0 + log(1d0-x)**2*log(x) 
     &+ log(1d0-x)*basis2_1(x) + 2*basis3_3(x) - 2*zeta3
      case(25)
         ris = -(pi**2*dlog(2d0))/12d0 + dlog(2d0)**3/6d0 
     &+ (pi**2*log(1d0-x))/12d0 - (dlog(2d0)**2*log(1d0-x))/2d0 
     &+ (dlog(2d0)*log(1d0-x)**2)/2d0 - basis3_6(x) + (7*zeta3)/8d0
      case(26)
         ris = (pi**2*log(1d0-x))/6d0 - basis3_3(x) + zeta3
      case(27)
         ris = -log(1d0-x)**3/6d0
         
      end  select


      endif
      HPL3=ris
      return
      end
      

c     #####################################################
c     basis3_1(z) = cli3(z) 
c     basis3_2(z) = cli3(-z)
c     basis3_3(z) = cli3(1-z)
c     basis3_4(z) = cli3(1/(1+z)) 
c     basis3_5(z) = cli3((1+z)/2) 
c     basis3_6(z) = cli3((1-z)/2) 
c     basis3_7(z) = cli3((1-z)/(1+z)) 
c     basis3_8(z) = cli3(2z/(z-1))
c     #####################################################

c ---------------------------------------------------------
c     basis3_1(z) = cli3(z) 
      double complex function basis3_1(x)
      implicit none
      double complex x,cli3           
      basis3_1=cli3(x)
      return
      end
c ---------------------------------------------------------
c     basis3_2(z) = cli3(-z)
      double complex function basis3_2(x)
      implicit none
      double complex x,cli3
      basis3_2=cli3(-x)
      return
      end
c ---------------------------------------------------------
c     basis3_3(z) = cli3(1-z)
      double complex function basis3_3(x)
      implicit none
      double complex x,cli3
      basis3_3 = cli3(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis3_4(z) = cli3(1/(1+z)) 
      double complex function basis3_4(x)
      implicit none
      double complex x,cli3
      basis3_4 = cli3(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_5(z) = cli3((1+z)/2) 
      double complex function basis3_5(x)
      implicit none
      double complex x,cli3
      basis3_5 = cli3((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_6(z) = cli3((1-z)/2) 
      double complex function basis3_6(x)
      implicit none
      double complex x,cli3
      basis3_6 = cli3((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_7(z) = cli3((1-z)/(1+z)) 
      double complex function basis3_7(x)
      implicit none
      double complex x,cli3
      basis3_7 = cli3((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_8(z) = cli3(2z/(z-1))
      double complex function basis3_8(x)
      implicit none
      double complex x,cli3
      basis3_8 = cli3(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------

C---  The single basis function we need to calculate: Li3
      
      double complex  function cli3(z)
      implicit none
      double complex ris, z, bsli3_inside,bsli3_outside, wcli3
      double precision zabs,border, pi, zeta2, zeta3
      
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
     
      border = 0.3d0
      zabs = abs(z)

      if (z.eq.(dcmplx(1d0,0d0))) then 
         ris = dcmplx(zeta3,0d0)
      else
         
         if (zabs.le.1d0) then
            if (zabs.le.border) then 
               ris=bsli3_inside(z)
            else
               ris=bsli3_outside(z)
            endif
         else
            ris=wcli3(1d0/z)-log(-z)**3/6d0-zeta2*log(-z)
         endif
      endif
      cli3=ris
      return
      end
      
      double complex  function wcli3(z)
      implicit none
      double complex z, cli3
      wcli3 =  cli3(z)
      return
      end
   
   
C---- expansion of trilogarithm in y = - log(1-z) with Bernoulli numbers  
C------requires  routine fbern3 in coefficients.F for the  coefficients 
C-------of the series  expansion 
      
      double  complex function bsli3_inside(z)
      implicit none
      integer i, Nmin
      double complex elem, ris, z, zb 
      double precision prec, fbern3, coeffi
      prec = 1d-8 
      Nmin=60
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmin 
         coeffi=fbern3(i)
         elem = zb**(i+1)*coeffi/(i+1)
         ris = ris+elem
      enddo
      
      bsli3_inside=ris 
      return 
      end

C---- expansion of the trilogarithm in log(z) with Zeta values  
C------requires  routine zetaval3 in coefficients.F for the  coefficients 
C-------of the series  expansion 
C-------- used for border < |z| < 1
      
      double  complex function bsli3_outside(z)
      implicit none
      integer i, Nmax
      double complex elem, ris, z, zb, coeffi
      double precision zetaval3
      Nmax=60
      zb = log(z)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmax 
         if (i.eq.2) then
            coeffi = (1d0 + 0.5d0 -log(-zb))/2d0
         else
            coeffi = dcmplx(zetaval3(i),0d0)
         endif
         elem = zb**i*coeffi
         ris = ris+elem
      enddo
      
      bsli3_outside=ris 
      return 
      end
