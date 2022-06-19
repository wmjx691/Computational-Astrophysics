! Test problem for 3D MHD
!
!      
! not ready     
      subroutine initial(den,px,py,pz,ene,bx,by,bz)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      real*8 pg,rhoe,rsq,B0,d0,v0,p0,pi
      real*8 az(ibeg-2:iend+2,jbeg-2:jend+2)
      pi = 4.d0*datan(1.d0)
c      B0 = 5.d0
      write(*,*)'warnnig!: this test has not setup the correct initial '
      write(*,*)'          condition.'
      do k=kbeg-kbuf, kend+kbuf
      do j=jbeg-jbuf, jend+jbuf
      do i=ibeg-ibuf, iend+ibuf

          den(i,j,k)= 1.d0
           px(i,j,k)= -dsin(2.d0*pi*y(j))       
           py(i,j,k)=  dsin(2.d0*pi*x(i))       
           pz(i,j,k)= 0.d0       
           bx(i,j,k)= -dsin(2.d0*pi*y(j))/gam
           by(i,j,k)=  dsin(4.d0*pi*x(i))/gam      
           bz(i,j,k)= 0.d0       
           pg     = 1.d0/gam   ! pressure
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
      enddo 
      enddo 
      enddo
      end
