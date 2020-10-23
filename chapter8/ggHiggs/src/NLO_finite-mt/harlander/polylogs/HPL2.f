C---- Evaluation of rank 2 HPLs with complex argument 

c --- Real part of HPL2     
      double precision function HPL2real(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2real = dreal(HPL2(n1,n2,x))
      return
      end

c --- Imaginary part of HPL2     
      double precision function HPL2im(n1,n2,xr,xi)
      implicit none
      double precision xr,xi
      integer n1,n2
      double complex x,HPL2
      x=dcmplx(xr,xi)
      HPL2im = dimag(HPL2(n1,n2,x))
      return
      end


C---  maps to basis function (dilogarithm)

      double  complex function HPL2(n1,n2, z)
      implicit none 
      double precision pi
      double complex z, ris,myi
c      double complex cli2 ! not needed ATM, due to mappings to basis functions
      double complex basis2_1,basis2_2,basis2_3
      integer n1,n2,j
      
      pi=3.1415926535897932385D0
      myi=dcmplx(0d0,1d0)
      ris=dcmplx(10d70,10d70)

      j = 3*(n1+1) + (n2+1) +1
      
      if (z.eq.dcmplx(0d0,0d0)) then
         if (j.ne.5) then
            ris = 0d0
         else
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL2: "
            print*, "HPL3(",n1,",",n2,",",z
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         endif
      else if (z.eq.dcmplx(1d0,0d0)) then
         if (j.eq.7.or.j.eq.9) then
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL2: "
            print*, "HPL3(",n1,",",n2,",",z
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop
         else
            select case(j)
            case(1)
               ris = dlog(2d0)**2/2d0
            case(2)
               ris = -pi**2/12d0
            case(3)
               ris = pi**2/12d0 - dlog(2d0)**2/2d0
            case(4)
               ris = pi**2/12d0
            case(5)
               ris = 0d0
            case(6)
               ris = pi**2/6d0
            case(8)
               ris = -pi**2/6d0
            end select
         endif
      else if (z.eq.dcmplx(-1d0,0d0)) then
         if (j.le.3) then
            print*, ""
            print*, "****************"
            print*, "ERROR in HPL2: "
            print*, "HPL3(",n1,",",n2,",",z
     &           ,") is divergent!"
            print*, "Aborting..."
            print*,"****************"
            stop   
         else
            select case(j)
            case(4)
               ris = -pi**2/6d0
            case(5)
               ris = -pi**2/2d0
            case(6)
               ris = -pi**2/12d0
            case(7)
               ris = pi**2/12d0 - dlog(2d0)**2/2d0
            case(8)
               ris = pi**2/12d0 - myi*pi*dlog(2d0)
            case(9)
               ris = dlog(2d0)**2/2d0
            end select
         endif
      else

c ########## hack to get the branch cuts right #################
         if (dimag(z).eq.0d0) then
            z = z + dcmplx(0d0,1d-60)
         endif
c ##############################################################

         select case(j)

c     #####################################################
c     basis2_1(z) = cli2(z)
c     basis2_2(z) = cli2(-z)) 
c     basis2_3(z) = cli2((1-z)/2)
c     #####################################################

      case(1)
         ris=log(1d0 + z)**2/2d0
      case(2)
         ris=basis2_2(z) + log(z)*log(1d0 + z)
      case(3)
         ris= pi**2/12d0 - dlog(2d0)**2/2d0 + dlog(2d0)*log(1d0-z) 
     &        - log(1d0-z)*log(1d0+z) - basis2_3(z)
      case(4)
         ris=-basis2_2(z)
      case(5)
         ris=log(z)**2/2d0
      case(6)
         ris=basis2_1(z)
      case(7)
         ris=-pi**2/12d0 + basis2_3(z) + dlog(2d0)**2/2d0 
     &        -dlog(2d0)*log(1d0 - z)
      case(8)
         ris=-basis2_1(z)-log(1d0 - z)*log(z)
      case(9)
         ris=log(1d0 - z)**2/2d0
      end select
      endif
      HPL2=ris 
      return
      end

c ---------------------------------------------------------
      double complex function basis2_1(x)
      implicit none
      double complex x,xcli2            
      basis2_1=xcli2(x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_2(x)
      implicit none
      double complex x,xcli2
      basis2_2=xcli2(-x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_3(x)
      implicit none
      double complex x,xcli2
      basis2_3=xcli2((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c ---------------------------------------------------------
c ---------------------------------------------------------


C--- maps to convergent region

      double complex  function xcli2(z)
      implicit none
      double complex ris, z, bsli2_inside,bsli2_outside, wcli2
      double precision zabs, pi, zeta2, border

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0

      border = 0.3d0 
      zabs = abs(z)
      if (z.eq.(dcmplx(1d0,0d0))) then 
         ris = dcmplx(zeta2,0d0)
      else
         if (zabs.le.1d0) then
            if (zabs.le.border) then 
               ris=bsli2_inside(z)
            else
               ris=bsli2_outside(z)
            endif
         else
            ris=-wcli2(1d0/z)-zeta2-0.5d0*log(-z)**2
         endif
      endif
         xcli2=ris
         return
      end
      
c---  fake recursion
      
      double complex  function wcli2(z)
      implicit none
      double complex z, xcli2
      wcli2 =  xcli2(z)
      return
      end


C---- expansion of dilogarithm in y = - log(1-z) with Bernoulli numbers  
C------ requires  routine fbern in bernoulli.F for the  coefficients 
C------- of the series  expansion 

      double  complex function bsli2_inside(z)
      implicit none
      integer i, Nmin
      double complex elem, ris, z, zb 
      double precision prec, fbern
      prec = 1d-8 
      Nmin=60
      zb = dcmplx(1d0,0d0)-z
      zb = -log(zb)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmin
      elem = zb**(i+1)*fbern(i)/(i+1)
      ris = ris+elem
      enddo
      bsli2_inside=ris 
      return 
      end

C---- expansion of the dilogarithm in log(z) with Zeta values  
C------requires  routine zetaval2 in coefficients.F for the  coefficients 
C-------of the series  expansion 
C-------- used for border < |z| < 1
      
      double  complex function bsli2_outside(z)
      implicit none
      integer i, Nmax
      double complex elem, ris, z, zb, coeffi
      double precision zetaval2
      Nmax=60
      zb = log(z)
      ris = dcmplx(0d0, 0d0)
      do i=0,Nmax 
         if (i.eq.1) then
            coeffi = 1d0 -log(-zb)
         else
            coeffi = dcmplx(zetaval2(i),0d0)
         endif
         elem = zb**i*coeffi
         ris = ris+elem
      enddo
      
      bsli2_outside=ris 
      return 
      end

C--------              ----------
C----- HPLs  of Rank 1  

      double  complex function HPL1(n1, z)
      implicit none
      double complex z, ris, x
      integer n1
      x=0d0
      ris = log(x)   
      if (abs(n1).gt.1) then
         print*, "Error in HPL1: Index out of range"
      else
         select case(n1)
      case(-1)
         ris =log(1.0d0 + z)
      case(0)
         ris=log(z)
      case(1)
         ris=-log(1.0d0 - z)
      end select
      endif
      HPL1=ris 
      return
      end
