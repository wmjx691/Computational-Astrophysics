      subroutine divclean(dt,dx,dy,dz,
     &                    den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                    den2,px2,py2,pz2,ene2,bx2,by2,bz2)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'divFree.h'
      real*8 monopole,maxmonopole,avemonopole
      real*8 btemp1,btemp2

      do k=kbeg,kend
       do j=jbeg,jend
        do i=ibeg,iend
         vxt2(i,j,k) = px2(i,j,k)/den2(i,j,k)
         vyt2(i,j,k) = py2(i,j,k)/den2(i,j,k)
c         vxt(i,j,k)  = px(i,j,k)/den(i,j,k)
c         vyt(i,j,k)  = py(i,j,k)/den(i,j,k)
          if (D3D)then
           vzt2(i,j,k) = pz2(i,j,k)/den2(i,j,k)
c           vzt(i,j,k)  = pz(i,j,k)/den(i,j,k)
          endif        
        enddo
       enddo
      enddo
      vxt = px/den
      vyt = py/den
      if (D3D)then
      vzt = pz/den
      endif
!---  for periodic b.c.
      call divcleanBdry(vxt,vyt,vzt,vxt2,vyt2,vzt2,bx2,by2,bz2)
!---  Keep ini  
c      call divcleanBdryini(vxt,vyt,vzt,vxt2,vyt2,vzt2,
c     &                     bx, by, bz, bx2, by2, bz2)
!---  Div clean (field-CD)  
!-------------------------------------------------------------
!                For 2D      
!-------------------------------------------------------------
      if (D2D) then
       k=1
       do j=jbeg-1,jend+1
        do i=ibeg-1,iend+1
         Ez(i,j,k) = -0.5d0*(vxt2(i,j,k)*by2(i,j,k) -
     &                       vyt2(i,j,k)*bx2(i,j,k) +
     &                       vxt(i,j,k) *by(i,j,k)  -
     &                       vyt(i,j,k) *bx(i,j,k)   )  
        enddo
       enddo
!---  Update new B field
       do j=jbeg,jend
        do i=ibeg,iend
         btemp(i,j,k) = bx2(i,j,k)**2+by2(i,j,k)**2+bz2(i,j,k)**2 !B*
         bx2(i,j,k)=bx(i,j,k)-0.5d0*(Ez(i,j+1,k)-Ez(i,j-1,k))*dt/dy
         by2(i,j,k)=by(i,j,k)+0.5d0*(Ez(i+1,j,k)-Ez(i-1,j,k))*dt/dx
        enddo
       enddo

!---  Update new total energy ?
       !do j=jbeg,jend
        !do i=ibeg,iend
        ! btemp1 = bx2(i,j,k)**2+by2(i,j,k)**2+bz2(i,j,k)**2
        ! ene2(i,j,k)=ene(i,j,k) + 0.5d0*(btemp1-btemp(i,j,k))
        !enddo
       !enddo


!---  Calculate the div B for testing code
      maxmonopole=0.d0
       do j=jbeg+1,jend-1
        do i=ibeg+1,iend-1
          monopole=   0.5d0*(-bx2(i-1,j,k)+bx2(i+1,j,k))/dx
     &              + 0.5d0*(-by2(i,j-1,k)+by2(i,j+1,k))/dy
          maxmonopole = dmax1(maxmonopole,monopole)
        enddo
       enddo
      write(*,*)'Max Div B =',maxmonopole

!-------------------------------------------------------------
!                For 3D      
!-------------------------------------------------------------
      elseif(D3D) then
!---  Div clean (field-CD)  
      do k=kbeg-1,kend+1
       do j=jbeg-1,jend+1
        do i=ibeg-1,iend+1
         Ex(i,j,k) = -0.5d0*(vyt2(i,j,k)*bz2(i,j,k) -
     &                       vzt2(i,j,k)*by2(i,j,k) +
     &                       vyt(i,j,k) *bz(i,j,k)  -
     &                       vzt(i,j,k) *by(i,j,k)   )  
         Ey(i,j,k) = -0.5d0*(vzt2(i,j,k)*bx2(i,j,k) -
     &                       vxt2(i,j,k)*bz2(i,j,k) +
     &                       vzt(i,j,k) *bx(i,j,k)  -
     &                       vxt(i,j,k) *bz(i,j,k)   )  
         Ez(i,j,k) = -0.5d0*(vxt2(i,j,k)*by2(i,j,k) -
     &                       vyt2(i,j,k)*bx2(i,j,k) +
     &                       vxt(i,j,k) *by(i,j,k)  -
     &                       vyt(i,j,k) *bx(i,j,k)   )  
        enddo
       enddo
      enddo
!---  Update new B field
      do k=kbeg,kend
       do j=jbeg,jend
        do i=ibeg,iend
         bx2(i,j,k)=bx(i,j,k)-0.5d0*(Ez(i,j+1,k)-Ez(i,j-1,k))*dt/dy
     &                       +0.5d0*(Ey(i,j,k+1)-Ey(i,j,k-1))*dt/dz
         by2(i,j,k)=by(i,j,k)+0.5d0*(Ez(i+1,j,k)-Ez(i-1,j,k))*dt/dx
     &                       -0.5d0*(Ex(i,j,k+1)-Ex(i,j,k-1))*dt/dz
         bz2(i,j,k)=bz(i,j,k)-0.5d0*(Ey(i+1,j,k)-Ey(i-1,j,k))*dt/dx
     &                       +0.5d0*(Ex(i,j+1,k)-Ex(i,j-1,k))*dt/dy
        enddo
       enddo
      enddo
!---  Calculate the div B for testing code
      maxmonopole=0.d0
      avemonopole=0.d0
c      do k=kbeg+1,kend-1
c       do j=jbeg+1,jend-1
c        do i=ibeg+1,iend-1
      do k=kbeg,kend
       do j=jbeg,jend
        do i=ibeg,iend
          monopole=   0.5d0*(-bx2(i-1,j,k)+bx2(i+1,j,k))/dx
     &              + 0.5d0*(-by2(i,j-1,k)+by2(i,j+1,k))/dy
     &              + 0.5d0*(-bz2(i,j,k-1)+bz2(i,j,k+1))/dz
          avemonopole = avemonopole + monopole
          maxmonopole = dmax1(maxmonopole,monopole)
        enddo
       enddo
      enddo
      write(*,*)'Max DivB=',maxmonopole
      write(*,*)'Ave DivB=',(avemonopole/(iend*jend*kend))
      endif
      end
