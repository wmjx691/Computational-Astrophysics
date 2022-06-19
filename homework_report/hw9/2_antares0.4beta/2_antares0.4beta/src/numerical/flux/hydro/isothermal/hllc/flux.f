C================================================================
C
C     HLLC Flux Solver isothermal
C
C     KCPAN ASIAA
C
C     Apr. 2007
C================================================================
      SUBROUTINE FLUX(UL,UR,slopeL,slopeR,FX,dt,sndL,sndR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'fluid.h'
      include 'parameters'

      DIMENSION  FX(5),FL(5),FR(5),FCs(5),Fs(5)
      DIMENSION  UL(5),UR(5),ULs(5),URs(5), FsL(5), FsR(5)
      dimension  slopeL(5), slopeR(5)
      real*8     spdL, spdR, sndL, sndR
      REAL*8     uxL,uyL,uzL
      REAL*8     uxR,uyR,uzR
      REAL*8     rhoL,rhoR,pTL,pTR
      REAL*8     usqr2,bsqL,bsqR
      REAL*8     tempL,tempR,asq,AA,Xtemp
      REAL*8     cf_L,cf_R,dlamdaf_nL,dlamdaf_pL
      REAL*8     dlamdaf_nR,dlamdaf_pR
      REAL*8     SL,SR,SLs,SRs
      REAL*8     rho_hll,px_hll
      REAL*8     Fhll_rho,Fhll_px
      REAL*8     rhos,uxs,pycs,pzcs,bycs,bzcs,bxc
      real*8     eneL, eneR
      real*8     UL, UR, SM

      sndL=iso_snd
      sndR=iso_snd

      AA=1.d0
      PT=0.d0

      rhoL = UL(1) + 0.5d0*slopeL(1)
      uxL  = UL(2) + 0.5d0*slopeL(2)
      uyL  = UL(3) + 0.5d0*slopeL(3)
      uzL  = UL(4) + 0.5d0*slopeL(4)
c      eneL = UL(5) + 0.5d0*slopeL(5)

      rhoR = UR(1) - 0.5d0*slopeR(1)
      uxR  = UR(2) - 0.5d0*slopeR(2)
      uyR  = UR(3) - 0.5d0*slopeR(3)
      uzR  = UR(4) - 0.5d0*slopeR(4)
c      eneR = UR(5) - 0.5d0*slopeR(5)

      pTL = sndL*sndL*rhoL
      pxL = uxL*rhoL
      pyL = uyL*rhoL
      pzL = uzL*rhoL

      if(pTL .le. 0.d0) then
         rhoL = UL(1)
         uxL  = UL(2)
         uyL  = UL(3)
         uzL  = UL(4)
c         eneL = UL(5)

         pTL = sndL*sndL*rhoL
         pxL = uxL*rhoL
         pyL = uyL*rhoL
         pzL = uzL*rhoL
c        write(*,*)  'ooooops L'
      endif

      pTR = sndR*sndR*rhoR
      pxR = uxR*rhoR
      pyR = uyR*rhoR
      pzR = uzR*rhoR

      if(pTR .le. 0.d0) then
        rhoR = UR(1)
        uxR  = UR(2)
        uyR  = UR(3)
        uzR  = UR(4)
c        eneR = UR(5)

        pTR = sndR*sndR*rhoR
        pxR = uxR*rhoR
        pyR = uyR*rhoR
        pzR = uzR*rhoR
c        write(*,*) 'oooooooops R'
      endif

       cf_L = sndL
       cf_R = sndR
       dlamdaf_pL = uxL + cf_L
       dlamdaf_nL = uxL - cf_L 
       dlamdaf_pR = uxR + cf_R
       dlamdaf_nR = uxR - cf_R 

       SL = dmin1(dlamdaf_nL,dlamdaf_nR)
       SR = dmax1(dlamdaf_pL,dlamdaf_pR)

       SM = ((SR-uxR)*rhoR*uxR-(SL-uxL)*rhoL*uxL-pTR+pTL)
     &      /((SR-uxR)*rhoR-(SL-uxL)*rhoL)

       pstar = rhoL*(uxL-SL)*(uxL-SM)+pTL

      rhosL = rhoL*(SL-uxL)/(SL-SM)
c       pxsL = ((SL-uxL)*rhoL*uxL+(pstar-pTL))/(SL-SM)
c       pysL = (SL-uxL)*rhoL*vxL/(SL-SM)
      pxsL = rhosL*SM
      pysL = rhosL*uyL
      pzsL = rhosL*uzL
c      enesL = ((SL-uxL)*eneL-pTL*uxL+pstar*SM)/(SL-SM)

      rhosR = rhoR*(SR-uxR)/(SR-SM)
c       pxsR = ((SR-uxR)*rhoR*uxR+(pstar-pTR))/(SR-SM)
c       pysR = (SR-uxR)*rhoR*vxR/(SR-SM)
      pxsR = rhosR*SM
      pysR = rhosR*uyR
      pzsR = rhosR*uzR
c      enesR = ((SR-uxR)*eneR-pTR*uxR+pstar*SM)/(SR-SM)

       FL(1) = rhoL*uxL
       FL(2) = rhoL*uxL**2.d0 +pTL 
       FL(3) = rhoL*uxL*uyL
       FL(4) = rhoL*uxL*uzL
c       FL(5) = (eneL+pTL)*uxL

       FR(1) = rhoR*uxR
       FR(2) = rhoR*uxR**2.d0 +pTR 
       FR(3) = rhoR*uxR*uyR 
       FR(4) = rhoR*uxR*uzR
c       FR(5) = (eneR+pTR)*uxR

       
       Fs(1)= (SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL+PT) 
       Fs(2)= (SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))
     &        /(SR-SL+PT) 
       Fs(3)= (SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))
     &        /(SR-SL+PT)
       Fs(4)= (SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))
     &        /(SR-SL+PT)
c       Fs(5)= (SR*FL(5)-SL*FR(5)+SR*SL*(UR(5)-UL(5)))/(SR-SL+PT) 

      FsL(1)= FL(1) + SL*(rhosL-rhoL)
      FsL(2)= FL(2) + SL*(pxsL-pxL)
      FsL(3)= FL(3) + SL*(pysL-pyL)
      FsL(4)= FL(4) + SL*(pzsL-pzL)
c      FsL(5)= FL(5) + SL*(enesL-eneL)

      FsR(1)= FR(1) + SR*(rhosR-rhoR)
      FsR(2)= FR(2) + SR*(pxsR-pxR)
      FsR(3)= FR(3) + SR*(pysR-pyR)
      FsR(4)= FR(4) + SR*(pzsR-pzR)
c      FsR(5)= FR(5) + SR*(enesR-eneR)


      IF (SL .gt. 0.d0)then
        do k = 1,4!5
        FX(k)=FL(k)
        enddo
      ELSEIF ( SL .le. 0.d0 .and. SM .ge. 0.d0)then
        do k = 1,4
        FX(k)=FsL(k)
        enddo
      ELSEIF ( SM .le. 0.d0 .and. SR .ge. 0.d0)then
        do k = 1,4
        FX(k)=FsR(k)
        enddo
      ELSE
        do k = 1,4
        FX(k)=FR(k)
        enddo
      ENDIF
      END

