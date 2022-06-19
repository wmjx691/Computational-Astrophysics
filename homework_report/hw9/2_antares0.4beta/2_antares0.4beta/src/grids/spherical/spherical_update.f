      subroutine spherical_update(dtdx,coef,cx,cy,cz,nd,S_den,
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
      real*8  rsq,rsqL,rsqR,ftempL,ftempR,denL,denR
 
      k=1
      j=1
c       do k=kbeg,kend
c       do j=jbeg,jend
       do i=ibeg,iend
       rsq  = x(i)**2 
       rsqL = xl(i-1)**2   ! x(i-1/2)
       rsqR = xl(i)**2     ! x(i+1/2)
       denL = 0.5d0*(den(i-1,j,k)+den(i,j,k))
       denR = 0.5d0*(den(i+1,j,k)+den(i,j,k))
       ftempL =F_px(i-cx,j-cy,k-cz)-(F_den(i-cx,j-cy,k-cz))**2/denL
       ftempR =F_px(i   ,j   ,k   )-(F_den(i   ,j   ,k   ))**2/denR

       den3(i,j,k)=den2(i,j,k)+coef(1)*
     &         (rsqL*F_den(i-cx,j-cy,k-cz)-rsqR*F_den(i,j,k))*dtdx/rsq
        px3(i,j,k)= px2(i,j,k)+coef(2)*
     &         (rsqL*(F_px(i-cx,j-cy,k-cz)-ftempL)-
     &          rsqR*(F_px(i,j,k)         -ftempR))*dtdx/rsq
     &         + coef(2)*(ftempL-ftempR)*dtdx

        py3(i,j,k)= py2(i,j,k)!+coef(3)*
c     &         (F_py(i-cx,j-cy,k-cz)-F_py(i,j,k))*dtdx

        pz3(i,j,k)= pz2(i,j,k)!+coef(4)*
c     &         (F_pz(i-cx,j-cy,k-cz)-F_pz(i,j,k))*dtdx

        if (.not. ISOTHERMAL)then
       ene3(i,j,k)=ene2(i,j,k)+coef(5)*
     &         (rsqL*F_ene(i-cx,j-cy,k-cz)-rsqR*F_ene(i,j,k))*dtdx/rsq
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
c       enddo
c       enddo

       end
