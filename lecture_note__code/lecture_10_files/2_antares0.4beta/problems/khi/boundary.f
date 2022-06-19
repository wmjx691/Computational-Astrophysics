      subroutine boundary(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                    den2,px2,py2,pz2,ene2,bx2,by2,bz2,uini)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      real*8 uini(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf,1:8)

!=================================================================
!                       Boundary Condition        
!=================================================================

!       inner :
! X
        call bnd_inner_periodic(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                          den2,px2,py2,pz2,ene2,bx2,by2,bz2,1)
c        call bnd_inner_reflect(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                     den2,px2,py2,pz2,ene2,bx2,by2,bz2,1,uini)
! Y
        if (D2D .or. D3D)then
        call bnd_inner_reflect(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                          den2,px2,py2,pz2,ene2,bx2,by2,bz2,2)
c        call bnd_inner_ini(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                     den2,px2,py2,pz2,ene2,bx2,by2,bz2,2,uini)
       endif
! Z 
       if (D3D)then
        call bnd_inner_outflow(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                          den2,px2,py2,pz2,ene2,bx2,by2,bz2,3)
c        call bnd_inner_ini(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                     den2,px2,py2,pz2,ene2,bx2,by2,bz2,3,uini)
       endif

!      outer :

! X       
        call bnd_outter_periodic(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                           den2,px2,py2,pz2,ene2,bx2,by2,bz2,1)
c        call bnd_outter_outflow(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,1,uini)
! Y
        if (D2D .or. D3D)then
        call bnd_outter_reflect(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                           den2,px2,py2,pz2,ene2,bx2,by2,bz2,2)
c        call bnd_outter_ini(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,2,uini)
       endif
! Z
       if (D3D)then
        call bnd_outter_outflow(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                           den2,px2,py2,pz2,ene2,bx2,by2,bz2,3)
c        call bnd_outter_ini(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
c     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,3,uini)
       endif


      end
