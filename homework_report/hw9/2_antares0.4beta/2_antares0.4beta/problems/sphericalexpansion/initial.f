! Test problem for 2D MHD
!
! Spherical Explosion      
      
      subroutine initial(den,px,py,pz,ene,bx,by,bz)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      real*8 pg,rhoe,rsq,B0,d0,v0,p0,pi
      real*8 az(ibeg-2:iend+2,jbeg-2:jend+2)
      pi = 4.d0*datan(1.d0)
      B0 = 5.d0
      
      do k=kbeg-kbuf, kend+kbuf
      do j=jbeg-jbuf, jend+jbuf
      do i=ibeg-ibuf, iend+ibuf
        if(D3D)then
          rsq = x(i)**2 +y(j)**2 +z(k)**2
        elseif(D2D)then
          rsq = x(i)**2 +y(j)**2
        endif
          r   = dsqrt(rsq)
        if ( r .le. 10.d0) then
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0
           py(i,j,k)= 0.d0
           pz(i,j,k)= 0.d0       
           bx(i,j,k)= 0.d0
           by(i,j,k)= B0/dsqrt(pi)
           bz(i,j,k)= 0.d0       
           pg     = 100.d0
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
        else      
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0       
           py(i,j,k)= 0.d0       
           pz(i,j,k)= 0.d0       
           bx(i,j,k)= 0.d0       
           by(i,j,k)= B0/dsqrt(pi)      
           bz(i,j,k)= 0.d0       
           pg     = 1.d0
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
        endif
      enddo 
      enddo 
      enddo
      end
