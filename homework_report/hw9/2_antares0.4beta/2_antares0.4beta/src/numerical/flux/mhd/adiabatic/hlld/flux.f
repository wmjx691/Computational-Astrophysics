!===================================================================
!HLLC solver: Hydrodynamic flux solver 
!Reference: Takahiro Miyoshi & kanya Kusano, JCP 208, p315-344, 2005
!Five fluxes will be calculated:
!     input: 
!          UL(:) => left state
!          UR(:) => right state
!          gam   => ratio of heat capacity
!     output: F(:)=> FLUX 
!          F(1) => density 
!          F(2) => x-momentum
!          F(3) => y-momentum
!          F(4) => z-momentum
!          F(5) => total energy
!          F(6) => Bx
!          F(7) => By
!          F(8) => Bz
!===================================================================
! Feb. Evora, Portugal. Hsiang-Hsu Wang
!===================================================================

      subroutine FLUX(UL,UR,slopeL,slopeR, F, dt,sndL,sndR)
      implicit none
      include 'fluid.h'
      include 'parameters'
      real*8 UL(8), UR(8), F(8), UssL(8), UssR(8),FssL(8),FssR(8)
      real*8 UsL(8), UsR(8), FL(8), FR(8), Fs(8), FsL(8), FsR(8)
      real*8 slopeL(8),slopeR(8)
      real*8 rhoL, rhoR  ! density
      real*8 xmomL, xmomR, ymomL, ymomR, zmomL, zmomR! momentum
      real*8 BxL, ByL, BzL, BxR, ByR, BzR,Bx ! magnetic field
      real*8 eneL, eneR  ! total energy
      real*8 pL, pR  ! left/right pressure
      real*8 eL, eR  ! left/right internal energy
      real*8 xuL, yuL, zuL, xuR, yuR, zuR  ! left/right fluid speed
      real*8 alfvenL, alfvenR  ! alfven's wave speed
      real*8 fastL, fastR, slowL, slowR ! fast/slow wave speed
      real*8 SL, SR  ! most negative/positive wave speed
      real*8 SM      ! speed of contact wave
      real*8 SLs, SRs ! SL*, SR*
      real*8 rhosL, rhosR !  rho*
      real*8 rhossL, rhossR !  rho**
      real*8 xusL,yusL,zusL, xusR,yusR,zusR ! u*, v*, w*
      real*8 xussL,yussL,zussL, xussR,yussR,zussR ! u**,v**,w**
      real*8 xmomsL, xmomsR !
      real*8 enesL, enesR !  e*
      real*8 enessL, enessR ! e**
      real*8 psT,pTL, pTR, pTsL, pTssL, pTsR, pTssR  ! psT = PT* = PT**
      real*8 BxsL, BysL, BzsL, BxsR, BysR, BzsR ! B*
      real*8 BxssL,ByssL,BzssL, BxssR, ByssR, BzssR ! B**
      real*8 BsqL, BsqR  ! B^2
      real*8 sndL,sndR
!===== make life easier    
       rhoL = UL(1)+ 0.5d0*slopeL(1)
        xuL = UL(2)+ 0.5d0*slopeL(2)
        yuL = UL(3)+ 0.5d0*slopeL(3)
        zuL = UL(4)+ 0.5d0*slopeL(4)
       eneL = UL(5)+ 0.5d0*slopeL(5)
        BxL = UL(6)+ 0.5d0*slopeL(6)
        ByL = UL(7)+ 0.5d0*slopeL(7)
        BzL = UL(8)+ 0.5d0*slopeL(8)

       rhoR = UR(1)- 0.5d0*slopeR(1)
        xuR = UR(2)- 0.5d0*slopeR(2)
        yuR = UR(3)- 0.5d0*slopeR(3)
        zuR = UR(4)- 0.5d0*slopeR(4)
       eneR = UR(5)- 0.5d0*slopeR(5)
        BxR = UR(6)- 0.5d0*slopeR(6)
        ByR = UR(7)- 0.5d0*slopeR(7)
        BzR = UR(8)- 0.5d0*slopeR(8)

       Bx = 0.5d0*(UL(6)+UR(6)) 
       !Bx = 0.5d0*(BxL+BxR) 
       !BxL = Bx
       !BxR = Bx

       BsqL = Bx**2+ByL**2+BzL**2
         eL = eneL-0.5d0*(rhoL*(xuL**2+yuL**2+zuL**2)+BsqL)    
         pL = (gam-1.d0)*eL
        pTL = pL + 0.5d0*BsqL
       alfvenL = Bx/dsqrt(rhoL)
       fastL = dsqrt((gam*pL+BsqL+dsqrt((gam*pL+BsqL)**2 
     &                 -4.d0*gam*pL*Bx**2))/(2.d0*rhoL))
       slowL = dsqrt((gam*pL+BsqL-dsqrt((gam*pL+BsqL)**2 
     &                 -4.d0*gam*pL*Bx**2))/(2.d0*rhoL))
 
       !write(*,*) "fastL, slowL", fastL, slowL

       BsqR = Bx**2+ByR**2+BzR**2
         eR = eneR-0.5d0*(rhoR*(xuR**2+yuR**2+zuR**2)+BsqR)
         pR = (gam-1.d0)*eR
        pTR = pR + 0.5d0*BsqR
       alfvenR = Bx/dsqrt(rhoR)
       fastR = dsqrt((gam*pR+BsqR+dsqrt((gam*pR+BsqR)**2 
     &                 -4.d0*gam*pR*Bx**2))/(2*rhoR))
       slowR = dsqrt((gam*pR+BsqR-dsqrt((gam*pR+BsqR)**2 
     &                 -4.d0*gam*pR*Bx**2))/(2*rhoR))

!====== HLLD calculation ==========
      SL = dmin1(xuL,xuR)-dmax1(fastL,fastR)
      SR = dmax1(xuL,xuR)+dmax1(fastL,fastR)
      SM = ((SR-xuR)*rhoR*xuR-(SL-xuL)*rhoL*xuL-pTR+pTL)/ 
     &      ((SR-xuR)*rhoR-(SL-xuL)*rhoL)
      psT = pTL + rhoL*(SL-xuL)*(SM-xuL)
!============== calculate star states 
      xusL  = SM
      xussL = SM
      xusR  = SM
      xussR = SM
      
      pTsL  = psT
      pTssL = psT
      pTsR  = psT
      pTssR = psT
! =========== left star state
      rhosL = rhoL*(SL-xuL)/(SL-SM)
      if(rhoL*(SL-xuL)*(SL-SM)-Bx**2 .eq. 0.d0) then
c      print *,"OoooopL!!!"
      yusL = yuL
      zusL = zuL
      BysL = ByL
      Bzsl = BzL      
      else
      yusL  = yuL-Bx*ByL*(SM-xuL)/(rhoL*(SL-xuL)*(SL-SM)-Bx**2)
      zusL  = zuL-Bx*BzL*(SM-xuL)/(rhoL*(SL-xuL)*(SL-SM)-Bx**2)
      BysL  = ByL*(rhoL*(SL-xuL)**2-Bx**2)/(rhoL*(SL-xuL)*(SL-SM)-Bx**2)
      BzsL  = BzL*(rhoL*(SL-xuL)**2-Bx**2)/(rhoL*(SL-xuL)*(SL-SM)-Bx**2)
      endif
      enesL = ((SL-xuL)*eneL-pTL*xuL+psT*SM+Bx*(xuL*Bx+yuL*ByL+zuL*BzL
     &          -xusL*Bx-yusL*BysL-zusL*BzsL))/(SL-SM)
! ========== right star state ================
      rhosR = rhoR*(SR-xuR)/(SR-SM)
      if(rhoR*(SR-xuR)*(SR-SM)-Bx**2 .eq. 0.d0) then
c      print *, "OoooopR!!!"
      yusR = yuR
      zusR = zuR
      BysR = ByR
      BzsR = BzR
      else
      yusR  = yuR-Bx*ByR*(SM-xuR)/(rhoR*(SR-xuR)*(SR-SM)-Bx**2)
      zusR  = zuR-Bx*BzR*(SM-xuR)/(rhoR*(SR-xuR)*(SR-SM)-Bx**2)
      BysR  = ByR*(rhoR*(SR-xuR)**2-Bx**2)/(rhoR*(SR-xuR)*(SR-SM)-Bx**2)
      BzsR  = BzR*(rhoR*(SR-xuR)**2-Bx**2)/(rhoR*(SR-xuR)*(SR-SM)-Bx**2)
      endif
      enesR = ((SR-xuR)*eneR-pTR*xuR+psT*SM+Bx*(xuR*Bx+yuR*ByR+zuR*BzR
     &           -xusR*Bx-yusR*BysR-zusR*BzsR))/(SR-SM)
!======= calculate StarStar state =======
      rhossL = rhosL
      rhossR = rhosR
         SLs = SM-dabs(Bx)/dsqrt(rhosL)
         SRs = SM+dabs(Bx)/dsqrt(rhosR)
     
!======= left starstar states ===========
      yussL=(dsqrt(rhosL)*yusL+dsqrt(rhosR)*yusR 
     &       +(BysR-BysL)*dsign(1.d0,Bx))/(dsqrt(rhosL)+dsqrt(rhosR))
      zussL=(dsqrt(rhosL)*zusL+dsqrt(rhosR)*zusR 
     &       +(BzsR-BzsL)*dsign(1.d0,Bx))/(dsqrt(rhosL)+dsqrt(rhosR))
      ByssL=(dsqrt(rhosL)*BysR+dsqrt(rhosR)*BysL+dsqrt(rhosL*rhosR) 
     &       *(yusR-yusL)*dsign(1.d0,Bx))/(dsqrt(rhosL)+dsqrt(rhosR))
      BzssL=(dsqrt(rhosL)*BzsR+dsqrt(rhosR)*BzsL+dsqrt(rhosL*rhosR) 
     &       *(zusR-zusL)*dsign(1.d0,Bx))/(dsqrt(rhosL)+dsqrt(rhosR))

      enessL=enesL-dsqrt(rhosL)*(xusL*Bx+yusL*BysL+zusL*BzsL   
     &       -xussL*Bx-yussL*ByssL-zussL*BzssL )*dsign(1.d0,Bx)
!======= right starstar states
      yussR=yussL
      zussR=zussL
      ByssR=ByssL
      BzssR=BzssL

      enessR=enesR+dsqrt(rhosR)*(xusR*Bx+yusR*BysR+zusR*BzsR   
     &       -xussR*Bx-yussR*ByssR-zussR*BzssR )*dsign(1.d0,Bx)
!======= compute flux =========
      UsL(1) = rhosL
      UsL(2) = rhosL*xusL
      UsL(3) = rhosL*yusL
      UsL(4) = rhosL*zusL
      UsL(5) = enesL
      UsL(6) = BxsL
      UsL(7) = BysL
      UsL(8) = BzsL

      UsR(1) = rhosR
      UsR(2) = rhosR*xusR
      UsR(3) = rhosR*yusR
      UsR(4) = rhosR*zusR
      UsR(5) = enesR
      UsR(6) = BxsR
      UsR(7) = BysR
      UsR(8) = BzsR

      UssL(1) = rhossL
      UssL(2) = rhossL*xussL
      UssL(3) = rhossL*yussL
      UssL(4) = rhossL*zussL
      UssL(5) = enessL
      UssL(6) = BxssL
      UssL(7) = ByssL
      UssL(8) = BzssL

      UssR(1) = rhossR
      UssR(2) = rhossR*xussR
      UssR(3) = rhossR*yussR
      UssR(4) = rhossR*zussR
      UssR(5) = enessR
      UssR(6) = BxssR
      UssR(7) = ByssR
      UssR(8) = BzssR

       FL(1) = rhoL*xuL
       FL(2) = rhoL*xuL**2+pTL-Bx**2
       FL(3) = rhoL*yuL*xuL-Bx*ByL
       FL(4) = rhoL*zuL*xuL-Bx*BzL
       FL(5) = (eneL+pTL)*xuL-Bx*(xuL*Bx+yuL*ByL+zuL*BzL)
       FL(6) = 0.d0
       FL(7) = ByL*xuL-Bx*yuL
       FL(8) = BzL*xuL-Bx*zuL

       FR(1) = rhoR*xuR
       FR(2) = rhoR*xuR**2+pTR-Bx**2
       FR(3) = rhoR*yuR*xuR-Bx*ByR
       FR(4) = rhoR*zuR*xuR-Bx*BzR
       FR(5) = (eneR+pTR)*xuR-Bx*(xuR*Bx+yuR*ByR+zuR*BzR)
       FR(6) = 0.d0
       FR(7) = ByR*xuR-Bx*yuR
       FR(8) = BzR*xuR-Bx*zuR

       FsL(1) = rhosL*xusL
       FsL(2) = rhosL*xusL**2+pTsL-Bx**2
       FsL(3) = rhosL*yusL*xusL-Bx*BysL
       FsL(4) = rhosL*zusL*xusL-Bx*BzsL
       FsL(5) = (enesL+pTsL)*xusL-Bx*(xusL*Bx+yusL*BysL+zusL*BzsL)     
       FsL(6) = 0.d0
       FsL(7) = BysL*xusL-Bx*yusL
       FsL(8) = BzsL*xusL-Bx*zusL
!       FsL = FL + SL*(UsL-UL)
!       FsR = FR + SR*(UsR-UR)
!       do i=1, 8
!          FsL(i) = FL(i) + SL*(UsL(i)-UL(i))
!          FsR(i) = FR(i) + SR*(UsR(i)-UR(i))  
!       enddo
       FsR(1) = rhosR*xusR
       FsR(2) = rhosR*xusR**2+pTsR-Bx**2 
       FsR(3) = rhosR*yusR*xusR-Bx*BysR
       FsR(4) = rhosR*zusR*xusR-Bx*BzsR
       FsR(5) = (enesR+pTsR)*xusR-Bx*(xusR*Bx+yusR*BysR+zusR*BzsR)
       FsR(6) = 0.d0
       FsR(7) = BysR*xusR-Bx*yusR
       FsR(8) = BzsR*xusR-Bx*zusR

!       FssL= FsL + SLs*(UssL-UsL)
!       FssR= FsR + SRs*(UssR-UsR) 
       FssL(1) = rhossL*xussL
       FssL(2) = rhossL*xussL**2+pTssL-Bx**2
       FssL(3) = rhossL*yussL*xussL-Bx*ByssL
       FssL(4) = rhossL*zussL*xussL-Bx*BzssL
       FssL(5) = (enessL+pTssL)*xussL-
     &            Bx*(xussL*Bx+yussL*ByssL+zussL*BzssL)
       FssL(6) = 0.d0
       FssL(7) = ByssL*xussL-Bx*yussL
       FssL(8) = BzssL*xussL-Bx*zussL

       FssR(1) = rhossR*xussR
       FssR(2) = rhossR*xussR**2+pTssR-Bx**2
       FssR(3) = rhossR*yussR*xussR-Bx*ByssR
       FssR(4) = rhossR*zussR*xussR-Bx*BzssR
       FssR(5) = (enessR+pTssR)*xussR-
     &            Bx*(xussR*Bx+yussR*ByssR+zussR*BzssR)
       FssR(6) = 0.d0
       FssR(7) = ByssR*xussR-Bx*yussR
       FssR(8) = BzssR*xussR-Bx*zussR
!====  choose correct flux
      IF (SL .gt. 0.d0)then
       do i=1,8
        F(i)=FL(i)
       enddo
      elseif( SL .le. 0.d0 .and. SLs .ge. 0.d0) then 
       do i=1,8
        F(i)=FsL(i)
       enddo
      elseif( SLs .le. 0.d0 .and. SM .ge. 0.d0) then
       do i=1,8
        F(i)=FssL(i)
       enddo
      elseif( SM .le. 0.d0 .and. SRs .ge. 0.d0) then
       do i=1,8
        F(i)=FssR(i)
       enddo
      elseif( SRs .le. 0.d0 .and. SR .ge. 0.d0) then
       do i=1,8
        F(i)=FsR(i)
       enddo
      else
       do i=1,8
        F(i)=FR(i)
       enddo
      ENDIF

      end
