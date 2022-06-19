! exact Riemann solver (isothermal)
!
! 2D isothermal
!      
      subroutine flux(VL,VR,slopeL,slopeR,Fx,dt,sndL,sndR)
      implicit double precision(a-h,o-z)
      include 'fluid.h'
      include 'parameters'
      dimension VL(5)
      dimension VR(5)
      dimension VM(5)
      dimension Fx(5)
      dimension slopeL(5)
      dimension slopeR(5)
      integer   n
      real*8    sndR,sndL,VL,VR,VM,slopeL,slopeR,Fx,c

      !write(*,*)'warnning:  exact Riemann solver has not modified '
      !write(*,*)'           for this version of antares codes'
      !write(*,*)'           try different flux solver.'
      !stop

      if(D3D)then
           write(*,*)'Error: this version of isothermal exact Riemaann'
           write(*,*)'       solver only for 1D or 2D cases'
           stop
      endif

      VL(1) = VL(1) + 0.5d0*slopeL(1)
      VL(2) = VL(2) + 0.5d0*slopeL(2)
      VL(3) = VL(3) + 0.5d0*slopeL(3)

      VR(1) = VR(1) - 0.5d0*slopeR(1)
      VR(2) = VR(2) - 0.5d0*slopeR(2)
      VR(3) = VR(3) - 0.5d0*slopeR(3)
      
      sndL=iso_snd
      sndR=iso_snd
      c   =iso_snd

      INtr = 10

      if (VL(1).eq.VR(1).and.VL(2).eq.VR(2)) then
           VM(1) = VR(1)
           VM(2) = VR(2)
       if(VL(2).ge. 0.   ) then
           VM(3) = VL(3)
       else
           VM(3) = VR(3)
       endif
       goto 1000
      endif
      rhoL = VL(1)
      rhoR = VR(1)
       uxL = VL(2)
       uxR = VR(2)
       uyL = VL(3)
       uyR = VR(3)

      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if( rhoL .ge. rhoR ) then
        if ((uxR-uxL) .gt. rarR) then
          goto 100
        elseif ((uxR-uxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else
        if ((uxR-uxL) .lt. shkR ) then
          goto 400
        elseif ((uxR-uxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
      endif

cccccc 1-rarefaction wave, 2-rarefaction wave
cccccc 1-rarefaction wave, 2-rarefaction wave
cccccc 1-rarefaction wave, 2-rarefaction wave
100   ustar   = 0.5d0*(uxL+uxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(ustar-uxL)/c)

      if ((uxL-c) .ge. 0.) then
         VM(2) = uxL
         VM(1) = rhoL
      elseif ((ustar-c) .ge. 0.) then
         VM(2) =  c
         VM(1) =  rhoL*dexp((uxL-c)/c)
      elseif ((ustar+c) .ge. 0.) then
         VM(2) =  ustar
         VM(1) =  rhostar
      elseif ((  uxR+c) .ge. 0.) then
         VM(2) =  -c
         VM(1) =  rhoR*dexp(-(uxR+c)/c)
      else
         VM(2) =   uxR
         VM(1) =  rhoR
      endif
      if ( VM(2) .ge. 0.) then
         VM(3) =  VL(3)
      else
         VM(3) =  VR(3)
      endif
      goto 1000
cccccc1-rarefaction wave, 2-shock wave
cccccc1-rarefaction wave, 2-shock wave
cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 -
     &  (uxR-uxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL))
     &  /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs(rho0-rhostar) .lt. 1.e-12) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rshostar,rhoL
      write(6,*) 'uxL    , uxR=',uxL,      uxR
      stop
210   ustar  = uxL - c*dlog(rhostar/rhoL)
      sigma2 = uxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         VM(2) =  uxR
         VM(1) = rhoR
      elseif ((ustar-c) .le. 0.) then
         VM(2) = ustar
         VM(1) = rhostar
      elseif ((uxL-c) .le. 0.) then
         VM(2) =  c
         VM(1) =  rhoL*dexp((uxL-c)/c)
      else
         VM(2) = uxL
         VM(1) = rhoL
      endif
      if ( VM(2) .ge. 0.) then
         VM(3) = VL(3)
      else
         VM(3) = VR(3)
      endif
      goto 1000
cccccc 1-shock wave, 2-rarefaction wave
cccccc 1-shock wave, 2-rarefaction wave
cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0-
     & (uxR-uxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR))
     &  /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs(rho0-rhostar) .lt. 1.e-12) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rshostar,rhoL
      write(6,*) 'uxL    , uxR=',uxL,      uxR
      stop
310   ustar  = uxR + c*dlog(rhostar/rhoR)
      sigma1 = uxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         VM(2) = uxL
         VM(1) = rhoL
      elseif ((ustar+c) .ge. 0.) then
         VM(2) = ustar
         VM(1) = rhostar
      elseif ((uxR+c) .ge. 0.) then
         VM(2) =  -c
         VM(1) =  rhoR*dexp(-(uxR+c)/c)
      else
         VM(2) = uxR
         VM(1) = rhoR
      endif
      if ( VM(2) .ge. 0.) then
         VM(3) = VL(3)
      else
         VM(3) = VR(3)
      endif
      goto 1000  
c       1-shock wave, 2-shock wave
c       1-shock wave, 2-shock wave
c       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = uxR-uxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      ustar    = uxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = uxL - c*dsqrt(rhostar/rhoL)
      sigma2   = uxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       VM(2) = uxL
       VM(1) = rhoL
      elseif (sigma2 .ge. 0.) then
       VM(2) = ustar
       VM(1) = rhostar
      else
       VM(2) = uxR
       VM(1) = rhoR
      endif
         if  (VM(2) .ge. 0.) then
       VM(3) = VL(3)
      else
       VM(3) = VR(3)
      endif
      goto 1000

1000  Fx(1) = VM(1)*(   VM(2)          )
      Fx(2) = VM(1)*( VM(2)*VM(2) + c*c)
      Fx(3) = VM(1)*( VM(2)*VM(3)      ) 
      Fx(4) = 0.d0
      Fx(5) = 0.d0 
      return
      end
