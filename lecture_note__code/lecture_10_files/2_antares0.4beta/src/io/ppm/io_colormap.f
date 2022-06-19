      
      subroutine colormap(color2)
      implicit none
      include 'color.h'
      integer*8 i,j,k,m

      open(3,file='color.palette')       
      read(3,*) color2
      end
      
