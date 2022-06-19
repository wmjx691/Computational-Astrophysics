      subroutine ddt(den,px,py,pz,ene,bx,by,bz,dx,dy,dt,dK)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      real*8 usq, uabs, dlamda, speed, vmax, tbd
      real*8 bsq,rhoe,prs,pstar,c,calf

      vmax=0.d0
      do k=kbeg,kend
      do j=jbeg,jend
      do i=ibeg,iend
       rho = den(i,j,k)
        vx = px(i,j,k)/rho
        vy = py(i,j,k)/rho
        vz = pz(i,j,k)/rho
        en = ene(i,j,k)
       usq = (vx*vx + vy*vy + vz*vz)
       uabs= dsqrt(usq) 
       if (MHD) then
       bsq = (bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2)
       else
       bsq = 0.d0
       endif       
       rhoe= en -0.5d0*(rho*usq +bsq)
       prs = rhoe*(gam-1.d0)
       pstar= prs + 0.5d0*bsq     
       dlamda= 0.5d0*rho/pstar
         if (dlamda .le. 0.d0 ) then
                 write(*,*)'negative lamda',i,j,k,dlamda
         endif
c        if (prs .le. 0.d0) then
c                write(*,*)'negative pressure',i,j,k,prs
c        endif        
       if (ISOTHERMAL)then
       c    = iso_snd
       else       
       c    = 1.d0/dsqrt(dlamda)*dsqrt(gam/2.d0)
c        c   = dsqrt(gam*prs/rho)       
       endif
       if (MHD)then
       calf = dsqrt(bsq/rho)
       else
       calf =0.d0
       endif       
c       speed = uabs + c + calf
       speed = dsqrt(uabs**2 + c**2 + calf**2)
       vmax  = max(Vmax,speed)
      enddo
      enddo
      enddo
      tbd= 1.d0/(vmax+1.d-10)
      dt = cfl*dmin1(10.d0,tbd)*dmin1(dx,dy)
      if (dt .le. smallt) then
              dt = smallt
      endif
      end
