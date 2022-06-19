      subroutine godnov3d(t,dt,dx,nd,
     &                    den ,px ,py ,pz, ene ,bx ,by ,bz , 
     &                    den2,px2,py2,pz2,ene2,bx2,by2,bz2, 
     &                    den3,px3,py3,pz3,ene3,bx3,by3,bz3) 
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'flux.h'
      include 'divFree.h'
      real*8 fslop
      external function fslop

      dtdx = dt/dx
      nx   = iend-ibeg+1
      ny   = jend-jbeg+1
      nz   = kend-kbeg+1

*--- nd ==1 -> x ,nd==2 -> y, nd==3 -> z
      if (nd .eq.1)then
        cx = 1
        cy = 0
        cz = 0

        f1 = 1
        f2 = 2
        f3 = 3
        f4 = 4
        f5 = 5
        f6 = 6
        f7 = 7
        f8 = 8

        coef(1)=1.d0
        coef(2)=1.d0
        coef(3)=1.d0
        coef(4)=1.d0
        coef(5)=1.d0
        coef(6)=1.d0
        coef(7)=1.d0
        coef(8)=1.d0

        g1 =1.d0
        g2 =1.d0
        g3 =1.d0
        g4 =1.d0
        g5 =1.d0
        g6 =1.d0
        g7 =1.d0
        g8 =1.d0

        do k=kbeg-cz,kend+cz
        do j=jbeg-cy,jend+cy
        do i=ibeg-cx,iend+cx
         rhoL=g1*den(i-cx,j-cy,k-cz)
         rhoM=g1*den(i   ,j   ,k   )
         rhoR=g1*den(i+cx,j+cy,k+cz)

          vxL=g2*px(i-cx,j-cy,k-cz)/rhoL
          vxM=g2*px(i   ,j   ,k   )/rhoM
          vxR=g2*px(i+cx,j+cy,k+cz)/rhoR

          vyL=g3*py(i-cx,j-cy,k-cz)/rhoL
          vyM=g3*py(i   ,j   ,k   )/rhoM
          vyR=g3*py(i+cx,j+cy,k+cz)/rhoR

          vzL=g4*pz(i-cx,j-cy,k-cz)/rhoL
          vzM=g4*pz(i   ,j   ,k   )/rhoM
          vzR=g4*pz(i+cx,j+cy,k+cz)/rhoR
          if (.not. ISOTHERMAL)then
          enL=g5*ene(i-cx,j-cy,k-cz)
          enM=g5*ene(i   ,j   ,k   )
          enR=g5*ene(i+cx,j+cy,k+cz)
          endif
          if (MHD) then
          bxL=g6*bx(i-cx,j-cy,k-cz)
          bxM=g6*bx(i   ,j   ,k   )
          bxR=g6*bx(i+cx,j+cy,k+cz)

          byL=g7*by(i-cx,j-cy,k-cz)
          byM=g7*by(i   ,j   ,k   )
          byR=g7*by(i+cx,j+cy,k+cz)

          bzL=g8*bz(i-cx,j-cy,k-cz)
          bzM=g8*bz(i   ,j   ,k   )
          bzR=g8*bz(i+cx,j+cy,k+cz)
          endif
          S_den(i,j,k)= fslop(rhoL,rhoM,rhoR)
          S_vx(i,j,k) = fslop( vxL, vxM, vxR)
          S_vy(i,j,k) = fslop( vyL, vyM, vyR)
          S_vz(i,j,k) = fslop( vzL, vzM, vzR)
          if (.not. ISOTHERMAL)then
          S_ene(i,j,k)= fslop( enL, enM, enR)
          endif
          if(MHD) then
          S_bx(i,j,k) = fslop( bxL, bxM, bxR)
          S_by(i,j,k) = fslop( byL, byM, byR)
          S_bz(i,j,k) = fslop( bzL, bzM, bzR)
          endif
        enddo  
        enddo
        enddo
        do k=jbeg-cz,kend
        do j=jbeg-cy,jend
        do i=ibeg-cx,iend
         UL(1)=g1*den(i,j,k)
         UL(2)=g2*px(i,j,k)
         UL(3)=g3*py(i,j,k)
         UL(4)=g4*pz(i,j,k)
         VL(1)=UL(1)
         VL(2)=UL(2)/UL(1)
         VL(3)=UL(3)/UL(1)
         VL(4)=UL(4)/UL(1)
         UR(1)=g1*den(i+cx,j+cy,k+cz)
         UR(2)=g2*px(i+cx,j+cy,k+cz)
         UR(3)=g3*py(i+cx,j+cy,k+cz)
         UR(4)=g4*pz(i+cx,j+cy,k+cz)
         VR(1)=UR(1)
         VR(2)=UR(2)/UR(1)
         VR(3)=UR(3)/UR(1)
         VR(4)=UR(4)/UR(1)
         slopeL(1)= S_den(i,j,k)
         slopeL(2)= S_vx(i,j,k)
         slopeL(3)= S_vy(i,j,k)
         slopeL(4)= S_vz(i,j,k)
         slopeR(1)= S_den(i+cx,j+cy,k+cz)
         slopeR(2)= S_vx(i+cx,j+cy,k+cz)
         slopeR(3)= S_vy(i+cx,j+cy,k+cz)
         slopeR(4)= S_vz(i+cx,j+cy,k+cz)
         if (.not. ISOTHERMAL)then
         UL(5)=g5*ene(i,j,k)
         VL(5)=UL(5)
         UR(5)=g5*ene(i+cx,j+cy,k+cz)
         VR(5)=UR(5)
         slopeL(5)= S_ene(i,j,k)
         slopeR(5)= S_ene(i+cx,j+cy,k+cz)
         endif
         if (MHD) then
         UL(6)=g6*bx(i,j,k)
         UL(7)=g7*by(i,j,k)
         UL(8)=g8*bz(i,j,k)
         VL(6)=UL(6)
         VL(7)=UL(7)
         VL(8)=UL(8)
         UR(6)=g6*bx(i+cx,j+cy,k+cz)
         UR(7)=g7*by(i+cx,j+cy,k+cz)
         UR(8)=g8*bz(i+cx,j+cy,k+cz)
         VR(6)=UR(6)
         VR(7)=UR(7)
         VR(8)=UR(8)
         slopeL(6)= S_bx(i,j,k)
         slopeL(7)= S_by(i,j,k)
         slopeL(8)= S_bz(i,j,k)
         slopeR(6)= S_bx(i+cx,j+cy,k+cz)
         slopeR(7)= S_by(i+cx,j+cy,k+cz)
         slopeR(8)= S_bz(i+cx,j+cy,k+cz)
         endif


         call flux(VL,VR,slopeL,slopeR,F,dt,snd,snd)

         F_den(i,j,k)=F(f1)
         F_px(i,j,k) =F(f2)
         F_py(i,j,k) =F(f3)
         F_pz(i,j,k) =F(f4)
         if (.not. ISOTHERMAL) then
         F_ene(i,j,k)=F(f5)
         endif
         if (MHD) then
         F_bx(i,j,k) =F(f6)
         F_by(i,j,k) =F(f7)
         F_bz(i,j,k) =F(f8)
         endif
        enddo
        enddo
        enddo
      elseif(nd .eq. 2)then
        cx = 0
        cy = 1
        cz = 0
        f1 = 1
        f2 = 3
        f3 = 2
        f4 = 4
        f5 = 5
        f6 = 7
        f7 = 6
        f8 = 8
        coef(1)= 1.d0
        coef(2)=-1.d0
        coef(3)= 1.d0
        coef(4)= 1.d0
        coef(5)= 1.d0
        coef(6)=-1.d0
        coef(7)= 1.d0
        coef(8)= 1.d0
        g1 = 1.d0
        g2 = 1.d0
        g3 =-1.d0
        g4 = 1.d0
        g5 = 1.d0
        g6 = 1.d0
        g7 =-1.d0
        g8 = 1.d0
        do k=kbeg-cz,kend+cz
        do j=jbeg-cy,jend+cy
        do i=ibeg-cx,iend+cx
         rhoL=g1*den(i-cx,j-cy,k-cz)
         rhoM=g1*den(i   ,j   ,k   )
         rhoR=g1*den(i+cx,j+cy,k+cz)
! X Y change          
          vxL=g2*py(i-cx,j-cy,k-cz)/rhoL
          vxM=g2*py(i   ,j   ,k   )/rhoM
          vxR=g2*py(i+cx,j+cy,k+cz)/rhoR
! X Y change          
          vyL=g3*px(i-cx,j-cy,k-cz)/rhoL
          vyM=g3*px(i   ,j   ,k   )/rhoM
          vyR=g3*px(i+cx,j+cy,k+cz)/rhoR
          vzL=g4*pz(i-cx,j-cy,k-cz)/rhoL
          vzM=g4*pz(i   ,j   ,k   )/rhoM
          vzR=g4*pz(i+cx,j+cy,k+cz)/rhoR
          S_den(i,j,k)= fslop(rhoL,rhoM,rhoR)
          S_vx(i,j,k) = fslop( vxL, vxM, vxR)
          S_vy(i,j,k) = fslop( vyL, vyM, vyR)
          S_vz(i,j,k) = fslop( vzL, vzM, vzR)

          if (.not. ISOTHERMAL) then
          enL=g5*ene(i-cx,j-cy,k-cz)
          enM=g5*ene(i   ,j   ,k   )
          enR=g5*ene(i+cx,j+cy,k+cz)
          S_ene(i,j,k)= fslop( enL, enM, enR)
          endif
          if (MHD) then
! X Y change          
          bxL=g6*by(i-cx,j-cy,k-cz)
          bxM=g6*by(i   ,j   ,k   )
          bxR=g6*by(i+cx,j+cy,k+cz)
! X Y change          
          byL=g7*bx(i-cx,j-cy,k-cz)
          byM=g7*bx(i   ,j   ,k   )
          byR=g7*bx(i+cx,j+cy,k+cz)
          bzL=g8*bz(i-cx,j-cy,k-cz)
          bzM=g8*bz(i   ,j   ,k   )
          bzR=g8*bz(i+cx,j+cy,k+cz)
          S_bx(i,j,k) = fslop( bxL, bxM, bxR)
          S_by(i,j,k) = fslop( byL, byM, byR)
          S_bz(i,j,k) = fslop( bzL, bzM, bzR)
          endif
        enddo
        enddo  
        enddo  
        do k=kbeg-cz,kend
        do j=jbeg-cy,jend
        do i=ibeg-cx,iend
         UL(1)=g1*den(i,j,k)
! X Y change          
         UL(2)=g2*py(i,j,k)
         UL(3)=g3*px(i,j,k)
         UL(4)=g4*pz(i,j,k)
         UR(1)=g1*den(i+cx,j+cy,k+cz)
         UR(2)=g2*py(i+cx,j+cy,k+cz)
         UR(3)=g3*px(i+cx,j+cy,k+cz)
         UR(4)=g4*pz(i+cx,j+cy,k+cz)
         VL(1)=UL(1)
         VL(2)=UL(2)/UL(1)
         VL(3)=UL(3)/UL(1)
         VL(4)=UL(4)/UL(1)
         VR(1)=UR(1)
         VR(2)=UR(2)/UR(1)
         VR(3)=UR(3)/UR(1)
         VR(4)=UR(4)/UR(1)
         slopeL(1)= S_den(i,j,k)
         slopeL(2)= S_vx(i,j,k)
         slopeL(3)= S_vy(i,j,k)
         slopeL(4)= S_vz(i,j,k)
         slopeR(1)= S_den(i+cx,j+cy,k+cz)
         slopeR(2)= S_vx(i+cx,j+cy,k+cz)
         slopeR(3)= S_vy(i+cx,j+cy,k+cz)
         slopeR(4)= S_vz(i+cx,j+cy,k+cz)

         if (.not.ISOTHERMAL)then
         UL(5)=g5*ene(i,j,k)
         VL(5)=UL(5)
         UR(5)=g5*ene(i+cx,j+cy,k+cz)
         VR(5)=UR(5)
         slopeL(5)= S_ene(i,j,k)
         slopeR(5)= S_ene(i+cx,j+cy,k+cz)
         endif
         if(MHD) then
         UL(6)=g6*by(i,j,k)
         UL(7)=g7*bx(i,j,k)
         UL(8)=g8*bz(i,j,k)
         VL(6)=UL(6)
         VL(7)=UL(7)
         VL(8)=UL(8)
         UR(6)=g6*by(i+cx,j+cy,k+cz)
         UR(7)=g7*bx(i+cx,j+cy,k+cz)
         UR(8)=g8*bz(i+cx,j+cy,k+cz)
         VR(6)=UR(6)
         VR(7)=UR(7)
         VR(8)=UR(8)
         slopeL(6)= S_bx(i,j,k)
         slopeL(7)= S_by(i,j,k)
         slopeL(8)= S_bz(i,j,k)
         slopeR(6)= S_bx(i+cx,j+cy,k+cz)
         slopeR(7)= S_by(i+cx,j+cy,k+cz)
         slopeR(8)= S_bz(i+cx,j+cy,k+cz)
         endif



         call flux(VL,VR,slopeL,slopeR,F,dt,snd,snd)

         F_den(i,j,k)=F(f1)
         F_px(i,j,k) =F(f2)
         F_py(i,j,k) =F(f3)
         F_pz(i,j,k) =F(f4)
         if (.not.ISOTHERMAL)then
         F_ene(i,j,k)=F(f5)
         endif
         if (MHD)then
         F_bx(i,j,k) =F(f6)
         F_by(i,j,k) =F(f7)
         F_bz(i,j,k) =F(f8)
         endif

        enddo
        enddo
        enddo
      elseif(nd .eq. 3)then
        cx = 0
        cy = 0
        cz = 1

        f1 = 1
        f2 = 4
        f3 = 3
        f4 = 2
        f5 = 5
        f6 = 8
        f7 = 7
        f8 = 6

        coef(1)= 1.d0
        coef(2)=-1.d0
        coef(3)= 1.d0
        coef(4)= 1.d0
        coef(5)= 1.d0
        coef(6)=-1.d0
        coef(7)= 1.d0
        coef(8)= 1.d0
        g1 = 1.d0
        g2 = 1.d0
        g3 = 1.d0
        g4 =-1.d0
        g5 = 1.d0
        g6 = 1.d0
        g7 = 1.d0
        g8 =-1.d0
        do k=kbeg-cz,kend+cz
        do j=jbeg-cy,jend+cy
        do i=ibeg-cx,iend+cx
         rhoL=g1*den(i-cx,j-cy,k-cz)
         rhoM=g1*den(i   ,j   ,k   )
         rhoR=g1*den(i+cx,j+cy,k+cz)

! X Y change          
          vxL=g2*pz(i-cx,j-cy,k-cz)/rhoL
          vxM=g2*pz(i   ,j   ,k   )/rhoM
          vxR=g2*pz(i+cx,j+cy,k+cz)/rhoR
! X Y change          
          vyL=g3*py(i-cx,j-cy,k-cz)/rhoL
          vyM=g3*py(i   ,j   ,k   )/rhoM
          vyR=g3*py(i+cx,j+cy,k+cz)/rhoR
! X Y change          
          vzL=g4*px(i-cx,j-cy,k-cz)/rhoL
          vzM=g4*px(i   ,j   ,k   )/rhoM
          vzR=g4*px(i+cx,j+cy,k+cz)/rhoR
          S_den(i,j,k)= fslop(rhoL,rhoM,rhoR)
          S_vx(i,j,k) = fslop( vxL, vxM, vxR)
          S_vy(i,j,k) = fslop( vyL, vyM, vyR)
          S_vz(i,j,k) = fslop( vzL, vzM, vzR)

          if(.not. ISOTHERMAL)then
          enL=g5*ene(i-cx,j-cy,k-cz)
          enM=g5*ene(i   ,j   ,k   )
          enR=g5*ene(i+cx,j+cy,k+cz)
          S_ene(i,j,k)= fslop( enL, enM, enR)
          endif
          if (MHD) then
! X Y change          
          bxL=g6*bz(i-cx,j-cy,k-cz)
          bxM=g6*bz(i   ,j   ,k   )
          bxR=g6*bz(i+cx,j+cy,k+cz)
! X Y change          
          byL=g7*by(i-cx,j-cy,k-cz)
          byM=g7*by(i   ,j   ,k   )
          byR=g7*by(i+cx,j+cy,k+cz)

          bzL=g8*bx(i-cx,j-cy,k-cz)
          bzM=g8*bx(i   ,j   ,k   )
          bzR=g8*bx(i+cx,j+cy,k+cz)
          S_bx(i,j,k) = fslop( bxL, bxM, bxR)
          S_by(i,j,k) = fslop( byL, byM, byR)
          S_bz(i,j,k) = fslop( bzL, bzM, bzR)
          endif
        enddo
        enddo  
        enddo  
        do k=kbeg-cz,kend
        do j=jbeg-cy,jend
        do i=ibeg-cx,iend
         UL(1)=g1*den(i,j,k)
! X Y change          
         UL(2)=g2*pz(i,j,k)
         UL(3)=g3*py(i,j,k)
         UL(4)=g4*px(i,j,k)
         UR(1)=g1*den(i+cx,j+cy,k+cz)
         UR(2)=g2*pz(i+cx,j+cy,k+cz)
         UR(3)=g3*py(i+cx,j+cy,k+cz)
         UR(4)=g4*px(i+cx,j+cy,k+cz)
         VL(1)=UL(1)
         VL(2)=UL(2)/UL(1)
         VL(3)=UL(3)/UL(1)
         VL(4)=UL(4)/UL(1)
         VR(1)=UR(1)
         VR(2)=UR(2)/UR(1)
         VR(3)=UR(3)/UR(1)
         VR(4)=UR(4)/UR(1)
         slopeL(1)= S_den(i,j,k)
         slopeL(2)= S_vx(i,j,k)
         slopeL(3)= S_vy(i,j,k)
         slopeL(4)= S_vz(i,j,k)
         slopeR(1)= S_den(i+cx,j+cy,k+cz)
         slopeR(2)= S_vx(i+cx,j+cy,k+cz)
         slopeR(3)= S_vy(i+cx,j+cy,k+cz)
         slopeR(4)= S_vz(i+cx,j+cy,k+cz)

         if (.not. ISOTHERMAL)then
         UL(5)=g5*ene(i,j,k)
         UR(5)=g5*ene(i+cx,j+cy,k+cz)
         VL(5)=UL(5)
         VR(5)=UR(5)
         slopeL(5)= S_ene(i,j,k)
         slopeR(5)= S_ene(i+cx,j+cy,k+cz)
         endif
         if (MHD)then
! X Y change          
         UL(6)=g6*bz(i,j,k)
         UL(7)=g7*by(i,j,k)
         UL(8)=g8*bx(i,j,k)
         UR(6)=g6*bz(i+cx,j+cy,k+cz)
         UR(7)=g7*by(i+cx,j+cy,k+cz)
         UR(8)=g8*bx(i+cx,j+cy,k+cz)
         VL(6)=UL(6)
         VL(7)=UL(7)
         VL(8)=UL(8)
         VR(6)=UR(6)
         VR(7)=UR(7)
         VR(8)=UR(8)
         slopeL(6)= S_bx(i,j,k)
         slopeL(7)= S_by(i,j,k)
         slopeL(8)= S_bz(i,j,k)
         slopeR(6)= S_bx(i+cx,j+cy,k+cz)
         slopeR(7)= S_by(i+cx,j+cy,k+cz)
         slopeR(8)= S_bz(i+cx,j+cy,k+cz)
         endif
          
         call flux(VL,VR,slopeL,slopeR,F,dt,snd,snd)

         F_den(i,j,k)=F(f1)
         F_px(i,j,k) =F(f2)
         F_py(i,j,k) =F(f3)
         F_pz(i,j,k) =F(f4)
         if (.not. ISOTHERMAL) then
         F_ene(i,j,k)=F(f5)
         endif
         if (MHD)then
         F_bx(i,j,k) =F(f6)
         F_by(i,j,k) =F(f7)
         F_bz(i,j,k) =F(f8)
         endif

        enddo
        enddo
        enddo

      else
       write(*,*)'no such direction, nd=',nd 
       stop
      endif  
     
!===   Choose Coordinate ====================================
       if (cartesian)then
       call cartesian_update(dtdx,coef,cx,cy,cz,nd,
     &             F_den,F_px,F_py,F_pz,F_ene,F_bx,F_by,F_bz , 
     &                      den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                      den3,px3,py3,pz3,ene3,bx3,by3,bz3)
       elseif(spherical)then 
       call spherical_update(dtdx,coef,cx,cy,cz,nd,S_den,
     &             F_den,F_px,F_py,F_pz,F_ene,F_bx,F_by,F_bz , 
     &                      den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                      den3,px3,py3,pz3,ene3,bx3,by3,bz3)
       elseif(cylindrical)then
       call cylindrical_update(dtdx,coef,cx,cy,cz,nd,
     &             F_den,F_px,F_py,F_pz,F_ene,F_bx,F_by,F_bz , 
     &                      den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                      den3,px3,py3,pz3,ene3,bx3,by3,bz3)
       else
          write(*,*)'error: No such coordinate'
          stop
       endif
!============================================================       
      end
