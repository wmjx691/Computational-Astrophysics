C234567891123456789212345678931234567894123456789512345678961234567897123
C========================================================================
C                 Gas-Kinetic Scheme                                                  
C                 Flux Splitting Method
c                 Ideal MHD
C                                                            Mar. 2007
C                                                            K.C. Pan
c                                                            ASIAA
C------------------------------------------------------------------------
      SUBROUTINE FLUX(w1,w2,slopeL,slopeR,FX,dt,sndL,sndR) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Include 'fluid.h'
      Include 'parameters'      
      DIMENSION  w1(8),w2(8),fnp(8),fnn(8),Ff(8),Fe(8),FX(8)
      DIMENSION  Une(8),Unp(8),Unn(8)
      real*8     slopeL(8),slopeR(8)
      INTEGER    nk
      REAL*8     ux1,uy1,uz1,bx1,by1,bz1
      REAL*8     ux2,uy2,uz2,bx2,by2,bz2
      REAL*8     rho0,ux0,uy0,uz0,bx0,by0,bz0,ene0
      REAL*8     rho1,rho2,ene1,ene2
      REAL*8     sinak,cosak
      REAL*8     unk1,utk1,bnk1,btk1,un0,ut0,bn0,bt0
      REAL*8     unk2,utk2,bnk2,byk2
      REAL*8     rhoe1,rhoe2,p1,p2,pstar1,pstar2mpstar0
      REAL*8     DLAM1,DLAM2
      REAL*8     u1,u2,b1,b2,usqr1,usqr2,bsqr1,bsqr2
      REAL*8     sndL,sndR,ETA
      ETA = 0.5d0
      pi  = 4.d0*datan(1.d0)
      nk = 1
      rho1=w1(1) + 0.5d0*slopeL(1)
       ux1=w1(2) + 0.5d0*slopeL(2)
       uy1=w1(3) + 0.5d0*slopeL(3)
       uz1=w1(4) + 0.5d0*slopeL(4)
      ene1=w1(5) + 0.5d0*slopeL(5)
       bx1=w1(6) + 0.5d0*slopeL(6)
       by1=w1(7) + 0.5d0*slopeL(7)
       bz1=w1(8) + 0.5d0*slopeL(8)

      rho2=w2(1) - 0.5d0*slopeR(1)
       ux2=w2(2) - 0.5d0*slopeR(2)
       uy2=w2(3) - 0.5d0*slopeR(3)
       uz2=w2(4) - 0.5d0*slopeR(4)
      ene2=w2(5) - 0.5d0*slopeR(5)
       bx2=w2(6) - 0.5d0*slopeR(6)
       by2=w2(7) - 0.5d0*slopeR(7)
       bz2=w2(8) - 0.5d0*slopeR(8)

      usqr1= ux1**2.0 + uy1**2.0 + uz1**2.0
      bsqr1= bx1**2.0 + by1**2.0 + bz1**2.0 
      u1= dsqrt(usqr1)
      b1= dsqrt(bsqr1)

       cosak = 1.d0
       sinak = 0.d0
    
      unk1 = ux1 
      utk1 = uy1
      bnk1 = bx1
      btk1 = by1
               
      rhoe1 = ene1-0.5*(rho1*usqr1+bsqr1)
      p1    = (GAM-1.0)*rhoe1
      pstar1= p1 + 0.5*(bnk1**2.0 + btk1**2.0 + bz1**2.0)
      DLAM1 = 0.5*rho1/pstar1
c      write(*,*)'re,p,p*,DL',rhoe1,p1,pstar1,DLAM1

      IF (DLAM1 .le. 0.0)THEN
      write(*,*)'LAMDA 1 negative'
      rho1=w1(1) 
       ux1=w1(2) 
       uy1=w1(3) 
       uz1=w1(4) 
      ene1=w1(5) 
       bx1=w1(6) 
       by1=w1(7) 
       bz1=w1(8) 

      rho2=w2(1) 
       ux2=w2(2) 
       uy2=w2(3) 
       uz2=w2(4) 
      ene2=w2(5) 
       bx2=w2(6) 
       by2=w2(7) 
       bz2=w2(8) 

      usqr1= ux1**2.0 + uy1**2.0 + uz1**2.0
      bsqr1= bx1**2.0 + by1**2.0 + bz1**2.0 
      u1= dsqrt(usqr1)
      b1= dsqrt(bsqr1)

       cosak = 1.d0
       sinak = 0.d0
    
      unk1 = ux1 
      utk1 = uy1
      bnk1 = bx1
      btk1 = by1
               
      rhoe1 = ene1-0.5*(rho1*usqr1+bsqr1)
      p1    = (GAM-1.0)*rhoe1
      pstar1= p1 + 0.5*(bnk1**2.0 + btk1**2.0 + bz1**2.0)
      DLAM1 = 0.5*rho1/pstar1
      ENDIF

c      write(*,*)'re,p,p*,DL',rhoe1,p1,pstar1,DLAM1
      
      V10p = 0.5*derfc(-dsqrt(DLAM1)*unk1)
c      V10n = 0.5*derfc( dsqrt(DLAM1)*unk1)
      V11p = unk1*V10p + 0.5*exp(-DLAM1*unk1**2.0)/dsqrt(DLAM1*pi)
c      V11n = unk1*V10n - 0.5*exp(-DLAM1*unk1**2.0)/dsqrt(DLAM1*pi)
   
      p01 = pstar1 - bnk1**2.0
      
      Fnp(1) = V11p* rho1     + V10p* 0.d0
      Fnp(2) = V11p* rho1*ux1 + V10p*(cosak*pstar1-bx1*bnk1)
      Fnp(3) = V11p* rho1*uy1 + V10p*(sinak*pstar1-by1*bnk1) 
      Fnp(4) = V11p* rho1*uz1 + V10p*(-bz1*bnk1)
      Fnp(5) = V11p* (ene1+0.5*p01) + V10p*(-bnk1*(btk1*utk1)
     &         +0.5*p01*unk1)
      Fnp(6) = V11p* (-sinak*btk1)  + V10p*( sinak*bnk1*utk1)
      Fnp(7) = V11p* cosak*btk1     + V10p*(-cosak*bnk1*utk1)
      Fnp(8) = V11p* bz1            + V10p*(-bnk1*uz1)

      Unp(1) = V11p* 0.d0 + V10p* rho1
      Unp(2) = V11p* rho1 + V10p* 0.d0
      Unp(3) = V11p* 0.d0 + V10p* rho1*utk1
      Unp(4) = V11p* 0.d0 + V10p* rho1*uz1
      Unp(5) = V11p* 0.5*rho1*unk1 + V10p*(ene1-0.5*rho1*unk1**2.0)
      Unp(6) = V11p* 0.d0 + V10p* bnk1
      Unp(7) = V11p* 0.d0 + V10p* btk1
      Unp(8) = V11p* 0.d0 + V10p* bz1
     
      unk2 = cosak*ux2 + sinak*uy2
      utk2 =-sinak*ux2 + cosak*uy2
      bnk2 = cosak*bx2 + sinak*by2
      btk2 =-sinak*bx2 + cosak*by2
      
      usqr2= ux2**2.0 + uy2**2.0 + uz2**2.0
      bsqr2= bx2**2.0 + by2**2.0 + bz2**2.0 
      u2= dsqrt(usqr2)
      b2= dsqrt(bsqr2)
      
      rhoe2 = ene2-0.5*(rho2*usqr2+bsqr2)
      p2    = (GAM-1.0)*rhoe2
      pstar2= p2 + 0.5*(bnk2**2.0 + btk2**2.0 + bz2**2.0)
      DLAM2 = 0.5*rho2/pstar2
      
      IF (DLAM1 .le. 0.0)THEN
      write(*,*)'LAMDA 2 negative'
      stop
      ENDIF

c      V20p = 0.5*derfc(-dsqrt(DLAM2)*unk2)
      V20n = 0.5*derfc( dsqrt(DLAM2)*unk2)
c      V21p = unk2*V20p + 0.5*exp(-DLAM2*unk2**2.0)/dsqrt(DLAM2*pi)
      V21n = unk2*V20n - 0.5*exp(-DLAM2*unk2**2.0)/dsqrt(DLAM2*pi)

      p02 = pstar2 - bnk2**2.0
      
      Fnn(1) = V21n* rho2     + V20n* 0.d0
      Fnn(2) = V21n* rho2*ux2 + V20n*(cosak*pstar2-bx2*bnk2)
      Fnn(3) = V21n* rho2*uy2 + V20n*(sinak*pstar2-by2*bnk2) 
      Fnn(4) = V21n* rho2*uz2 + V20n*(-bz2*bnk2)
      Fnn(5) = V21n* (ene2+0.5*p02) + V20n*(-bnk2*(btk2*utk2)
     &         +0.5*p02*unk2)
      Fnn(6) = V21n* (-sinak*btk2)  + V20n*( sinak*bnk2*utk2)
      Fnn(7) = V21n* cosak*btk2     + V20n*(-cosak*bnk2*utk2)
      Fnn(8) = V21n* bz2            + V20n*(-bnk2*uz2)
      
      Unn(1) = V21n* 0.d0 + V20n* rho2
      Unn(2) = V21n* rho2 + V20n* 0.d0
      Unn(3) = V21n* 0.d0 + V20n* rho2*utk2
      Unn(4) = V21n* 0.d0 + V20n* rho2*uz2
      Unn(5) = V21n* 0.5*rho2*unk2 + V20n*(ene2-0.5*rho2*unk2**2.0)
      Unn(6) = V21n* 0.d0 + V20n* bnk2
      Unn(7) = V21n* 0.d0 + V20n* btk2
      Unn(8) = V21n* 0.d0 + V20n* bz2
    
        do k = 1,8
        Ff(k)  = Fnp(k)+Fnn(k)
        Une(k) = Unp(k)+Unn(k)
        enddo

      rho0= Une(1)
      un0 = Une(2)/rho0 
      ut0 = Une(3)/rho0
      uz0 = Une(4)/rho0
      ene0= Une(5)
      bn0 = Une(6)
      bt0 = Une(7)
      bz0 = Une(8)

       ux0= un0
       uy0= ut0
       bx0= bn0
       by0= bt0
      
      usqr= ux0**2.0 + uy0**2.0 + uz0**2.0
      bsqr= bx0**2.0 + by0**2.0 + bz0**2.0
      rhoe0 = ene0-0.5*(rho0*usqr+bsqr)
      p0    = (GAM-1.0)*rhoe0
      pstar0= p0 + 0.5*(bn0**2.0 + bt0**2.0 + bz0**2.0)
      
      Fe(1) =  rho0*un0   
      Fe(2) =  rho0*ux0*un0-bx0*bn0+ cosak*pstar0
      Fe(3) =  rho0*uy0*un0-by0*bn0+ sinak*pstar0
      Fe(4) =  rho0*uz0*un0-bz0*bn0
      Fe(5) =  un0*(ene0+pstar0)-bn0*(un0*bn0+ut0*bt0+uz0*bz0)
      Fe(6) =  sinak*(bn0*ut0-bt0*un0)
      Fe(7) =  cosak*(bt0*un0-bn0*ut0)
      Fe(8) =  bz0*un0 - bn0*uz0
        
      do k = 1,8
      FX(k)= (1.d0-ETA)*Fe(k) + ETA*Ff(k)
      enddo
      
      END

C234567891123456789212345678931234567894123456789512345678961234567897123
