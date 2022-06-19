! Test problem for 1D hydro
!
! Woodward-Colella Test Case
!
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

        if ( x(i) .le. 0.1d0*xmax) then
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0
           py(i,j,k)= 0.d0
           pz(i,j,k)= 0.d0       
          ene(i,j,k)= 2500.d0
        elseif ( x(i) .le. 0.9d0*xmax) then
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0
           py(i,j,k)= 0.d0
           pz(i,j,k)= 0.d0       
          ene(i,j,k)= 0.025d0
        else      
          den(i,j,k)= 1.d0
           px(i,j,k)= 0.d0       
           py(i,j,k)= 0.d0       
           pz(i,j,k)= 0.d0       
          ene(i,j,k)= 250.d0
        endif
      enddo 
      enddo 
      enddo
      end
