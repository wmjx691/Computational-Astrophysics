      subroutine cylindrical_update(dtdx,coef,cx,cy,cz,nd,
     &             F_den,F_px,F_py,F_pz,F_ene,F_bx,F_by,F_bz, 
     &                      den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                      den3,px3,py3,pz3,ene3,bx3,by3,bz3)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'flux.h'
      include 'divFree.h'
 
!---  Note: in cylindrical coordinate ---------------
!
!            X --> r
!            Y --> z
!
!-----------------------------------------------------      
      if (nd .eq. 1)then
      k=1
c       do k=kbeg,kend
       do j=jbeg,jend
       do i=ibeg,iend

       den3(i,j,k)=den2(i,j,k)+coef(1)*
     &         (F_den(i-cx,j-cy,k-cz)-F_den(i,j,k))*dtdx
        px3(i,j,k)= px2(i,j,k)+coef(2)*
     &         (F_px(i-cx,j-cy,k-cz)-F_px(i,j,k))*dtdx
        py3(i,j,k)= py2(i,j,k)+coef(3)*
     &         (F_py(i-cx,j-cy,k-cz)-F_py(i,j,k))*dtdx
        pz3(i,j,k)= pz2(i,j,k)+coef(4)*
     &         (F_pz(i-cx,j-cy,k-cz)-F_pz(i,j,k))*dtdx
        if (.not. ISOTHERMAL)then
       ene3(i,j,k)=ene2(i,j,k)+coef(5)*
     &         (F_ene(i-cx,j-cy,k-cz)-F_ene(i,j,k))*dtdx
        endif
        if (MHD)then
        bx3(i,j,k)= bx2(i,j,k)+coef(6)*
     &         (F_bx(i-cx,j-cy,k-cz)-F_bx(i,j,k))*dtdx
        by3(i,j,k)= by2(i,j,k)+coef(7)*
     &         (F_by(i-cx,j-cy,k-cz)-F_by(i,j,k))*dtdx
        bz3(i,j,k)= bz2(i,j,k)+coef(8)*
     &         (F_bz(i-cx,j-cy,k-cz)-F_bz(i,j,k))*dtdx
        endif
       enddo
       enddo
c       enddo
     
       elseif (nd.eq.2)then

       k=1
c       do k=kbeg,kend
       do j=jbeg,jend
       do i=ibeg,iend

       den3(i,j,k)=den2(i,j,k)+coef(1)*
     &         (F_den(i-cx,j-cy,k-cz)-F_den(i,j,k))*dtdx
        px3(i,j,k)= px2(i,j,k)+coef(2)*
     &         (F_px(i-cx,j-cy,k-cz)-F_px(i,j,k))*dtdx
        py3(i,j,k)= py2(i,j,k)+coef(3)*
     &         (F_py(i-cx,j-cy,k-cz)-F_py(i,j,k))*dtdx
        pz3(i,j,k)= pz2(i,j,k)+coef(4)*
     &         (F_pz(i-cx,j-cy,k-cz)-F_pz(i,j,k))*dtdx
        if (.not. ISOTHERMAL)then
       ene3(i,j,k)=ene2(i,j,k)+coef(5)*
     &         (F_ene(i-cx,j-cy,k-cz)-F_ene(i,j,k))*dtdx
        endif
        if (MHD)then
        bx3(i,j,k)= bx2(i,j,k)+coef(6)*
     &         (F_bx(i-cx,j-cy,k-cz)-F_bx(i,j,k))*dtdx
        by3(i,j,k)= by2(i,j,k)+coef(7)*
     &         (F_by(i-cx,j-cy,k-cz)-F_by(i,j,k))*dtdx
        bz3(i,j,k)= bz2(i,j,k)+coef(8)*
     &         (F_bz(i-cx,j-cy,k-cz)-F_bz(i,j,k))*dtdx
        endif
       enddo
       enddo
c       enddo
       endif
       end
