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
      real*8 utemp(1:7),xsub,ysub,zsub,dsx,dsy,dsz
      real*8 utemp2(1:7)
      integer nsub,ii,jj,kk,ntotal
      pi = 4.d0*datan(1.d0)
c      B0 = 5.d0
      nsub = 6
      if (D3D) then
              write(*,*)'orsag tang problem is a 2D problem'
              stop
      endif
      k=1
      do j=jbeg-jbuf, jend+jbuf
      do i=ibeg-ibuf, iend+ibuf

          utemp2 = 0.d0
          ntotal = 0
        do ii = 1, nsub
        do jj = 1, nsub
          dsx = (x(i)-xl(i))*2.d0/dble(nsub)
          dsy = (y(i)-yl(i))*2.d0/dble(nsub)

          xsub = xl(i) + dble(ii)*dsx
          ysub = yl(j) + dble(jj)*dsy

           utemp(1) = 2.778d0
           utemp(2)= -2.778d0*dsin(ysub)       
           utemp(3)=  2.778d0*dsin(xsub)       
           utemp(4)= 0.d0       
           utemp(5)= -dsin(ysub)
           utemp(6)=  dsin(2.d0*xsub)      
           utemp(7)= 0.d0       
          
           utemp2 = utemp2 + utemp
           ntotal = ntotal + 1
        enddo
        enddo   

          den(i,j,k) = utemp2(1)/dble(ntotal)
           px(i,j,k) = utemp2(2)/dble(ntotal)
           py(i,j,k) = utemp2(3)/dble(ntotal)
           pz(i,j,k) = utemp2(4)/dble(ntotal)
           bx(i,j,k) = utemp2(5)/dble(ntotal)
           by(i,j,k) = utemp2(6)/dble(ntotal)
           bz(i,j,k) = utemp2(7)/dble(ntotal)
!          den(i,j,k)= 1.d0
!           px(i,j,k)= -dsin(2.d0*pi*y(j))       
!           py(i,j,k)=  dsin(2.d0*pi*x(i))       
!           pz(i,j,k)= 0.d0       
!           bx(i,j,k)= -dsin(2.d0*pi*y(j))/gam
!           by(i,j,k)=  dsin(4.d0*pi*x(i))/gam      
!           bz(i,j,k)= 0.d0       
           pg     = 1.667d0   ! pressure
          rhoe    = pg/(gam-1.d0)
          ene(i,j,k)= rhoe + 
     &                0.5d0*(bx(i,j,k)**2+by(i,j,k)**2+bz(i,j,k)**2 +
     &             (px(i,j,k)**2+py(i,j,k)**2+pz(i,j,k)**2)/den(i,j,k))
      enddo 
      enddo 
      end
