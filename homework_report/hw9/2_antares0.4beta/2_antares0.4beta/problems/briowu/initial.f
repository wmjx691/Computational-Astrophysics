! Test problem for 1D MHD
!
! Brio-Wu shock tube      
!      
      subroutine initial(den,px,py,pz,ene,bx,by,bz)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      real*8 pg,rhoe,rsq,B0,d0,v0,p0,pi
      real*8 az(ibeg-2:iend+2,jbeg-2:jend+2)
      
      do k=kbeg-kbuf, kend+kbuf
      do j=jbeg-jbuf, jend+jbuf
      do i=ibeg-ibuf, iend+ibuf

        if ( x(i) .le. 0.d0) then
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0
           py(i,j,k)= 0.d0
           pz(i,j,k)= 0.d0       
           bx(i,j,k)= 0.75d0
           by(i,j,k)= 1.d0
           bz(i,j,k)= 0.d0       
           pg     = 1.d0
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
        else      
          den(i,j,k)= 0.125d0
           px(i,j,k)= 0.d0       
           py(i,j,k)= 0.d0       
           pz(i,j,k)= 0.d0       
           bx(i,j,k)= 0.75d0       
           by(i,j,k)=-1.d0      
           bz(i,j,k)= 0.d0       
           pg     = 0.1d0
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
        endif
      enddo 
      enddo 
      enddo
      end
