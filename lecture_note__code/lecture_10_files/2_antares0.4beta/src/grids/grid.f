      subroutine grid 
      implicit none
      include 'fluid.h'
      include 'parameters'
      
      nx = iend - ibeg +1
      dx = (xmax-xmin)/(iend-ibeg+1)
      ny = jend - jbeg +1
      dy = (ymax-ymin)/(jend-jbeg+1)
      nz = kend - kbeg +1
      dz = (zmax-zmin)/(kend-kbeg+1)

      do i = ibeg-ibuf-1,iend+ibuf
         xl(i) = xmin+i*dx
      enddo
      do i = ibeg-ibuf,iend+ibuf
         x(i) = 0.5d0*(xl(i-1)+xl(i))  
      enddo
      do j = jbeg-jbuf-1,jend+jbuf
         yl(j) = ymin+j*dy
      enddo
      do j = jbeg-jbuf,jend+jbuf
         y(j) = 0.5d0*(yl(j-1)+yl(j))  
      enddo
      do k = kbeg-kbuf-1,kend+kbuf
         zl(k) = zmin+k*dz
      enddo
      do k = kbeg-kbuf,kend+kbuf
         z(k) = 0.5d0*(zl(k-1)+zl(k))  
      enddo
      end 
