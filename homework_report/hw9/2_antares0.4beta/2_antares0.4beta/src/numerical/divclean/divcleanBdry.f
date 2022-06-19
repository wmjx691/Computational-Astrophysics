      subroutine divcleanBdry(vxt,vyt,vzt,vxt2,vyt2,vzt2,bx,by,bz)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'divFree.h'

!     For periodic B.C. -----      
      if (D2D)then
      k = 1
      do j = jbeg,jend
      bx(ibeg-1,j,k) = bx(iend,j,k)
      by(ibeg-1,j,k) = by(iend,j,k)
      vxt(ibeg-1,j,k)= vxt(iend,j,k)
      vyt(ibeg-1,j,k)= vyt(iend,j,k)
      vxt2(ibeg-1,j,k)= vxt2(iend,j,k)
      vyt2(ibeg-1,j,k)= vyt2(iend,j,k)
      bx(iend+1,j,k) = bx(ibeg,j,k)
      by(iend+1,j,k) = by(ibeg,j,k)
      vxt(iend+1,j,k)= vxt(ibeg,j,k)
      vyt(iend+1,j,k)= vyt(ibeg,j,k)
      vxt2(iend+1,j,k)= vxt2(ibeg,j,k)
      vyt2(iend+1,j,k)= vyt2(ibeg,j,k)
      enddo
      do i = ibeg,iend
      bx(i,jbeg-1,k) = bx(i,jend,k)
      by(i,jbeg-1,k) = by(i,jend,k)
      vxt(i,jbeg-1,k)= vxt(i,jend,k)
      vyt(i,jbeg-1,k)= vyt(i,jend,k)
      vxt2(i,jbeg-1,k)= vxt2(i,jend,k)
      vyt2(i,jbeg-1,k)= vyt2(i,jend,k)
      bx(i,jend+1,k) = bx(i,jbeg,k)
      by(i,jend+1,k) = by(i,jbeg,k)
      vxt(i,jend+1,k)= vxt(i,jbeg,k)
      vyt(i,jend+1,k)= vyt(i,jbeg,k)
      vxt2(i,jend+1,k)= vxt2(i,jbeg,k)
      vyt2(i,jend+1,k)= vyt2(i,jbeg,k)
      enddo
!---------------------------------------------------      
      elseif(D3D)then
      do k=kbeg,kend
      do j=jbeg,jend
      bx(ibeg-1,j,k) = bx(iend,j,k)
      by(ibeg-1,j,k) = by(iend,j,k)
      bz(ibeg-1,j,k) = bz(iend,j,k)
      vxt(ibeg-1,j,k)= vxt(iend,j,k)
      vyt(ibeg-1,j,k)= vyt(iend,j,k)
      vzt(ibeg-1,j,k)= vzt(iend,j,k)
      vxt2(ibeg-1,j,k)= vxt2(iend,j,k)
      vyt2(ibeg-1,j,k)= vyt2(iend,j,k)
      vzt2(ibeg-1,j,k)= vzt2(iend,j,k)
      bx(iend+1,j,k) = bx(ibeg,j,k)
      by(iend+1,j,k) = by(ibeg,j,k)
      bz(iend+1,j,k) = bz(ibeg,j,k)
      vxt(iend+1,j,k)= vxt(ibeg,j,k)
      vyt(iend+1,j,k)= vyt(ibeg,j,k)
      vzt(iend+1,j,k)= vzt(ibeg,j,k)
      vxt2(iend+1,j,k)= vxt2(ibeg,j,k)
      vyt2(iend+1,j,k)= vyt2(ibeg,j,k)
      vzt2(iend+1,j,k)= vzt2(ibeg,j,k)
      enddo
      enddo

      do k=kbeg,kend
      do i=ibeg,iend
      bx(i,jbeg-1,k) = bx(i,jend,k)
      by(i,jbeg-1,k) = by(i,jend,k)
      bz(i,jbeg-1,k) = bz(i,jend,k)
      vxt(i,jbeg-1,k)= vxt(i,jend,k)
      vyt(i,jbeg-1,k)= vyt(i,jend,k)
      vzt(i,jbeg-1,k)= vzt(i,jend,k)
      vxt2(i,jbeg-1,k)= vxt2(i,jend,k)
      vyt2(i,jbeg-1,k)= vyt2(i,jend,k)
      vzt2(i,jbeg-1,k)= vzt2(i,jend,k)
      bx(i,jend+1,k) = bx(i,jbeg,k)
      by(i,jend+1,k) = by(i,jbeg,k)
      bz(i,jend+1,k) = bz(i,jbeg,k)
      vxt(i,jend+1,k)= vxt(i,jbeg,k)
      vyt(i,jend+1,k)= vyt(i,jbeg,k)
      vzt(i,jend+1,k)= vzt(i,jbeg,k)
      vxt2(i,jend+1,k)= vxt2(i,jbeg,k)
      vyt2(i,jend+1,k)= vyt2(i,jbeg,k)
      vzt2(i,jend+1,k)= vzt2(i,jbeg,k)
      enddo
      enddo

      do j=jbeg,jend
      do i=ibeg,iend
      bx(i,j,kbeg-1) = bx(i,j,kend)
      by(i,j,kbeg-1) = by(i,j,kend)
      bz(i,j,kbeg-1) = bz(i,j,kend)
      vxt(i,j,kbeg-1) = vxt(i,j,kend)
      vyt(i,j,kbeg-1) = vyt(i,j,kend)
      vzt(i,j,kbeg-1) = vzt(i,j,kend)
      vxt2(i,j,kbeg-1) = vxt2(i,j,kend)
      vyt2(i,j,kbeg-1) = vyt2(i,j,kend)
      vzt2(i,j,kbeg-1) = vzt2(i,j,kend)
      bx(i,j,kend+1) = bx(i,j,kbeg)
      by(i,j,kend+1) = by(i,j,kbeg)
      bz(i,j,kend+1) = bz(i,j,kbeg)
      vxt(i,j,kend+1) = vxt(i,j,kbeg)
      vyt(i,j,kend+1) = vyt(i,j,kbeg)
      vzt(i,j,kend+1) = vzt(i,j,kbeg)
      vxt2(i,j,kend+1) = vxt2(i,j,kbeg)
      vyt2(i,j,kend+1) = vyt2(i,j,kbeg)
      vzt2(i,j,kend+1) = vzt2(i,j,kbeg)
      enddo
      enddo
      else
       write(*,*)'No such Dim for Div clean at Bdry'
       stop
      endif       
      end



      subroutine divcleanBdryini(vxt,vyt,vzt,vxt2,vyt2,vzt2,
     &                           bx,by,bz,bx2,by2,bz2)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'divFree.h'

!     For keep ini B.C. -----      
      if (D2D)then
        write(*,*)'keep ini B.C for div clean only work for 3D'
        write(*,*)'please modify the B.C. in divclean.f'
        stop
!---------------------------------------------------      
      elseif(D3D)then
      do k=kbeg,kend
      do j=jbeg,jend
      bx2(ibeg-1,j,k) = bx(ibeg-1,j,k)
      by2(ibeg-1,j,k) = by(ibeg-1,j,k)
      bz2(ibeg-1,j,k) = bz(ibeg-1,j,k)
      vxt2(ibeg-1,j,k)= vxt(ibeg-1,j,k)
      vyt2(ibeg-1,j,k)= vyt(ibeg-1,j,k)
      vzt2(ibeg-1,j,k)= vzt(ibeg-1,j,k)
      bx2(iend+1,j,k) = bx(iend+1,j,k)
      by2(iend+1,j,k) = by(iend+1,j,k)
      bz2(iend+1,j,k) = bz(iend+1,j,k)
      vxt2(iend+1,j,k)= vxt(iend+1,j,k)
      vyt2(iend+1,j,k)= vyt(iend+1,j,k)
      vzt2(iend+1,j,k)= vzt(iend+1,j,k)
      enddo
      enddo

      do k=kbeg,kend
      do i=ibeg,iend
      bx2(i,jbeg-1,k) = bx(i,jbeg-1,k)
      by2(i,jbeg-1,k) = by(i,jbeg-1,k)
      bz2(i,jbeg-1,k) = bz(i,jbeg-1,k)
      vxt2(i,jbeg-1,k)= vxt(i,jbeg-1,k)
      vyt2(i,jbeg-1,k)= vyt(i,jbeg-1,k)
      vzt2(i,jbeg-1,k)= vzt(i,jbeg-1,k)
      bx2(i,jend+1,k) = bx(i,jend+1,k)
      by2(i,jend+1,k) = by(i,jend+1,k)
      bz2(i,jend+1,k) = bz(i,jend+1,k)
      vxt2(i,jend+1,k)= vxt(i,jend+1,k)
      vyt2(i,jend+1,k)= vyt(i,jend+1,k)
      vzt2(i,jend+1,k)= vzt(i,jend+1,k)
      enddo
      enddo

      do j=jbeg,jend
      do i=ibeg,iend
      bx2(i,j,kbeg-1) = bx(i,j,kbeg-1)
      by2(i,j,kbeg-1) = by(i,j,kbeg-1)
      bz2(i,j,kbeg-1) = bz(i,j,kbeg-1)
      vxt2(i,j,kbeg-1) = vxt(i,j,kbeg-1)
      vyt2(i,j,kbeg-1) = vyt(i,j,kbeg-1)
      vzt2(i,j,kbeg-1) = vzt(i,j,kbeg-1)
      bx2(i,j,kend+1) = bx(i,j,kend+1)
      by2(i,j,kend+1) = by(i,j,kend+1)
      bz2(i,j,kend+1) = bz(i,j,kend+1)
      vxt2(i,j,kend+1) = vxt(i,j,kend+1)
      vyt2(i,j,kend+1) = vyt(i,j,kend+1)
      vzt2(i,j,kend+1) = vzt(i,j,kend+1)
      enddo
      enddo
      else
       write(*,*)'No such Dim for Div clean at Bdry'
       stop
      endif       
      end
