* --------------------------------------------
* Sample FORTRAN 77 program to make PPM / BMP
*
*  1998-Apr-03 (long long ago)
*      :
*  2005-Jul-22 (added P3)
*      :
*  2006 Dec-22 (added BMP)
*
*  Image array rgb(3,*,*) is filled in subroutine mkbitmap()
*  and
*  saved in PPM or BMP format in subroutine pixout().
*
*                                   K. Hayashi
* --------------------------------------------
*
* 2007-Nov-27 (modified for CFD-MHD codes)
*                                   Odie K.C. Pan

       subroutine io_ppm(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,
     &                   IOC,color2)
       implicit none
       include 'fluid.h'
       include 'parameters'
       include 'conserve.h'
       include 'io.h'
       include 'color.h'
       integer   ihpixf,jvpixf,khpixf
       character*1 rgb(3,iend+2*ibuf,jend+2*jbuf)

       ihpixf = iend+2*ibuf
       jvpixf = jend+2*jbuf       
       khpixf = kend+2*kbuf

       call makemap(den,rgb,nprnt,ihpixf,jvpixf,khpixf,color2)
       call  pixout(rgb,nprnt,ihpixf,jvpixf,khpixf)
       if (D3D)then
       call makemapXZ(den,rgb,nprnt,ihpixf,jvpixf,khpixf,color2)
       call  pixoutXZ(rgb,nprnt,ihpixf,jvpixf,khpixf)
       endif
       end
*==============================================================
       subroutine makemap(den,rgb,nprnt,ihpixf,jvpixf,khpixf,color2)
       implicit none
       include 'fluid.h'
       include 'parameters'
       include 'color.h'
       integer   ihpixf,jvpixf,khpixf
       character*1 rgb(3,ihpixf,jvpixf)
       real*8 den(-1:ihpixf-2,-1:jvpixf-2,-1:khpixf-2)
       integer nrho(ihpixf,jvpixf)
       integer   m,ii,jj
       real*8    rhomin,rhomax 
       real*8    drho,rhorang 
       integer ilog
!       kend=(khpixf-4)
       ilog = 2   !  1 for linear scale
                  !  2 for log    scale
       if(D3D)then           
*=========================================
*---   Linear ouput 3D
       if (ilog .eq. 1) then
         rhomin = 1.d10
         rhomax = 0.d0
         do j = -1,jvpixf-2
         do i = -1,ihpixf-2
           if (rhomin .gt. den(i,j,kend/2)) then
             rhomin = den(i,j,kend/2)
           endif
           if (rhomax .lt. den(i,j,kend/2)) then
             rhomax = den(i,j,kend/2)
           endif
         enddo
         enddo
*---  Fixed rhomin and rhomax ---
c       rhomin =0.9d0
c       rhomax =2.1d0
*--------------------------------
        rhorang=rhomax-rhomin
        drho   =rhorang/255.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
       nrho(i+2,j+2)= ANINT( (den(i,j,kend/2)-rhomin)/drho )      
         if (rhomax.eq.rhomin)then
               nrho(i+2,j+2)=128
         endif      
         if (nrho(i+2,j+2) .gt. 256)then
               nrho(i+2,j+2)=256
         elseif(nrho(i+2,j+2).lt.1)then 
               nrho(i+2,j+2)=1
         endif    
       enddo
       enddo

       do k=1,3
       do j=1,jvpixf
       do i=1,ihpixf
        m = nrho(j,i)
       rgb(k,j,i)= char(color2(k,m))
       enddo
       enddo
       enddo
*========================================
*      Log scale output 3D

       elseif (ilog .eq. 2) then


       rhomin = 1.d10
       rhomax = 0.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
         if (rhomin .gt. dlog(den(i,j,kend/2))) then
           rhomin = dlog(den(i,j,kend/2))
         endif
         if (rhomax .lt. dlog(den(i,j,kend/2))) then
           rhomax = dlog(den(i,j,kend/2))
         endif
       enddo
       enddo
*---  Fixed rhomin and rhomax ---
c       rhomin = 1.d0
c       rhomax = 10.d0
*--------------------------------

       rhorang=rhomax-rhomin
       drho   =rhorang/255.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
        nrho(i+2,j+2)= ANINT( (dlog(den(i,j,kend/2))-rhomin)/drho )    
         if (rhomax.eq.rhomin)then
               nrho(i+2,j+2)=128
         endif      
         if (nrho(i+2,j+2) .gt. 256)then
               nrho(i+2,j+2)=256
         elseif(nrho(i+2,j+2).lt.1)then 
               nrho(i+2,j+2)=1
         endif    
       enddo
       enddo

         do k=1,3
         do j=1,jvpixf
         do i=1,ihpixf
           m = nrho(j,i)
           rgb(k,j,i)= char(color2(k,m))
         enddo
         enddo
         enddo
       endif

       elseif(D2D) then
*=========================================
*---   Linear ouput 2D 
                  
       if (ilog .eq. 1) then
         rhomin = 1.d10
         rhomax = 0.d0
         do j = -1,jvpixf-2
         do i = -1,ihpixf-2
           if (rhomin .gt. den(i,j,kbeg)) then
             rhomin = den(i,j,kbeg)
           endif
           if (rhomax .lt. den(i,j,kbeg)) then
             rhomax = den(i,j,kbeg)
           endif
         enddo
         enddo
*---  Fixed rhomin and rhomax ---
       rhomin =0.9d0
       rhomax =2.1d0
*--------------------------------
        rhorang=rhomax-rhomin
        drho   =rhorang/255.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
       nrho(i+2,j+2)= ANINT( (den(i,j,kbeg)-rhomin)/drho )      
         if (rhomax.eq.rhomin)then
               nrho(i+2,j+2)=128
         endif      
         if (nrho(i+2,j+2) .gt. 256)then
               nrho(i+2,j+2)=256
         elseif(nrho(i+2,j+2).lt.1)then 
               nrho(i+2,j+2)=1
         endif    
       enddo
       enddo

       do k=1,3
       do j=1,jvpixf
       do i=1,ihpixf
        m = nrho(j,i)
       rgb(k,j,i)= char(color2(k,m))
       enddo
       enddo
       enddo
*========================================
*      Log scale output  2D

       elseif (ilog .eq. 2) then


       rhomin = 1.d10
       rhomax = 0.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
         if (rhomin .gt. dlog(den(i,j,kbeg))) then
           rhomin = dlog(den(i,j,kbeg))
         endif
         if (rhomax .lt. dlog(den(i,j,kbeg))) then
           rhomax = dlog(den(i,j,kbeg))
         endif
       enddo
       enddo
!       write(*,*)'rhomin =',rhomin
!       write(*,*)'rhomax =',rhomax
*---  Fixed rhomin and rhomax ---
!       rhomin =-1.1d0
!       rhomax = 1.1d0
*--------------------------------

       rhorang=rhomax-rhomin
       drho   =rhorang/255.d0

       do j = -1,jvpixf-2
       do i = -1,ihpixf-2
        nrho(i+2,j+2)= ANINT( (dlog(den(i,j,kbeg))-rhomin)/drho )    
         if (rhomax.eq.rhomin)then
               nrho(i+2,j+2)=128
         endif      
         if (nrho(i+2,j+2) .gt. 256)then
               nrho(i+2,j+2)=256
         elseif(nrho(i+2,j+2).lt.1)then 
               nrho(i+2,j+2)=1
         endif    
       enddo
       enddo

         do k=1,3
         do j=1,jvpixf
         do i=1,ihpixf
           m = nrho(j,i)
           rgb(k,j,i)= char(color2(k,m))
         enddo
         enddo
         enddo
       endif

       else
       endif        
       end




* --------------------------------------------
*
* Notes
* o With a parameter ipixout set at 1, 2 or others,
*   this subroutine will generate PPM-P6(binary), PPM-P3(text),
*   or BMP(24bit depth without color table).
*
* o Some parts follow DEC-FORTRAN that had been defacto-standard long ago.
*   Some compilers today may not accept if "ipixout" is not 2.
*
* o g77 (ver. 3.3.3) works for all three choices.
* o Recent intel compiler (ver. 9 or so) works for all three choices.
*
* --------------------------------------------
*
       subroutine pixout(rgb,nframe,ihpixf,jvpixf,khpixf)
       implicit none
* interface arg.
       integer ihpixf, jvpixf, khpixf
       character*1 rgb(3,ihpixf,jvpixf)      ! RGB data array
       integer nframe
* local
       character*12 fnameout
       integer i, j, k
       integer idummy, icnt
       character*14 frmtstr
       character*54 headmsw
       character*4  byt4
       character*2  byt2
* choices
       integer ipixout
       parameter(ipixout = 1) ! 1 / 2 / other= PPM6, PPM3, BMP(24bit)

       if (ipixout .EQ. 1) then

* PPM P6

         write(fnameout,'(''rXY'',i4.4,''.ppm'')') nframe ! name of PPM file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing PPM (P6) file : ', fnameout
* header
         write(2,'(''P6'', 2(1x,i4),'' 255 '',$)')         ! some compiler may not accept this line.
     &     ihpixf, jvpixf
* image data
         idummy = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i8.8,''A,$)'')') idummy     ! make output "format"
         write(2,fmt=frmtstr)                              ! some compiler may not accept this line.
     &     (((rgb(k,i,j),k=1,3),i=1,ihpixf),j=jvpixf,1,-1) ! here, j (vertical address) runs from top to bottom.
         close(2)

       else if (ipixout .EQ. 2) then

* PPM P3 ! rather "safer" choice for many Fortran compiler(s).

         write(fnameout,'(''rXY'',i4.4,''.ppm'')') nframe ! name of PPM file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing PPM (P3) file : ', fnameout
* header
         write(2,'(A)') 'P3'
         write(2,'(2(1x,i4),'' 255 '')')  ihpixf, jvpixf
         icnt = 0
* image data
         do j = jvpixf, 1, -1                              ! here, j (vertical address) runs from top to bottom.
         do i = 1, ihpixf, 1
         do k = 1, 3
           idummy = ichar(rgb(k,i,j))
           icnt = icnt + 4
           if (icnt .LT. 60) then
             write(2,fmt='(1x,i3,$)') idummy               ! mind "$" is not standard.
           else
             write(2,fmt='(1x,i3)') idummy
             icnt = 0
           endif
         enddo
         enddo
         enddo
         write(2,'(A)') ' '
         close(2)

       else

* BMP (24bit depth)... this part works only when width is multiple of 4.

         idummy = mod(ihpixf, 4)
         if (idummy .NE. 0) then
           write(*,*) 'width must be multiple of 4'
           stop
         endif

         write(fnameout,'(''rXY'',i4.4,''.bmp'')') nframe ! name of BMP file
         open(unit=2,file=fnameout,status='unknown')
         write(*,*) 'Now writing BMP(24bit) file : ', fnameout
* header 1 (file header ; 1--14 byte)
         headmsw( 1: 2) = 'BM'             ! declaring this is BMP file
         idummy = 54 + ihpixf * jvpixf * 3 ! total file size = header + data
         call num2bit4(idummy,byt4)
         headmsw( 3: 6) = byt4(1:4)
         idummy = 0                        ! may be 0
         call num2bit2(idummy,byt2)
         headmsw( 7: 8) = byt2(1:2)
         idummy = 0                        ! may be 0
         call num2bit2(idummy,byt2)
         headmsw( 9:10) = byt2(1:2)
         idummy = 54                       ! must be 54 : total length of header
         call num2bit4(idummy,byt4)
         headmsw(11:14) = byt4(1:4)
* header 2 (bit-map header ; 13--54 byte)
         idummy = 40                       ! must be 40 : length of bit-map header
         call num2bit4(idummy,byt4)
         headmsw(15:18) = byt4(1:4)
         idummy = ihpixf                   ! width
         call num2bit4(idummy,byt4)
         headmsw(19:22) = byt4(1:4)
         idummy = jvpixf                   ! height
         call num2bit4(idummy,byt4)
         headmsw(23:26) = byt4(1:4)
         idummy = 1                        ! must be 1
         call num2bit2(idummy,byt2)
         headmsw(27:28) = byt2(1:2)
         idummy = 24                       ! must be 24 : color depth in bit.
         call num2bit2(idummy,byt2)
         headmsw(29:30) = byt2(1:2)
         idummy = 0                        ! may be 0 : compression method index
         call num2bit4(idummy,byt4)
         headmsw(31:34) = byt4(1:4)
         idummy = 0                        ! may be 0 : file size if compressed
         call num2bit4(idummy,byt4)
         headmsw(35:38) = byt4(1:4)
         idummy = 0                        ! arbit. : pixel per meter, horizontal
         call num2bit4(idummy,byt4)
         headmsw(39:42) = byt4(1:4)
         idummy = 0                        ! arbit. : pixel per meter, vertical
         call num2bit4(idummy,byt4)
         headmsw(43:46) = byt4(1:4)
         idummy = 0                        ! may be 0 here : num. of color used
         call num2bit4(idummy,byt4)
         headmsw(47:50) = byt4(1:4)
         idummy = 0                        ! may be 0 here : num. of important color
         call num2bit4(idummy,byt4)
         headmsw(51:54) = byt4(1:4)

* writing header part
         write(2,'(a54,$)') headmsw(1:54)
* image data
         idummy = ihpixf * jvpixf * 3
         write(frmtstr,'(''('',i8.8,''A,$)'')') idummy
         write(2,fmt=frmtstr)
     &     (((rgb(k,i,j),k=3,1,-1),i=1,ihpixf),j=1,jvpixf) ! writing in BGR order, not RGB.
         close(2)

       endif

       return
       end subroutine pixout

* --------------------------------------
* convert number to 8-bit characters
* --------------------------------------

       subroutine num2bit4(inum,byt4)
       implicit none
       integer inum
       character*4 byt4
       integer idummy1, idummy2
       idummy1 = inum
       idummy2 = idummy1 / 256**3
       byt4(4:4) = char(idummy2)
       idummy1 =-idummy2 * 256**3 +idummy1
       idummy2 = idummy1 / 256**2
       byt4(3:3) = char(idummy2)
       idummy1 =-idummy2 * 256**2 +idummy1
       idummy2 = idummy1 / 256
       byt4(2:2) = char(idummy2)
       idummy1 =-idummy2 * 256    +idummy1
       byt4(1:1) = char(idummy1)
       return
       end subroutine num2bit4

* ------

       subroutine num2bit2(inum,byt2)
       implicit none
       integer inum
       character*2 byt2
       integer idummy1, idummy2
       idummy1 = inum
       idummy2 = idummy1 / 256
       byt2(2:2) = char(idummy2)
       idummy1 =-idummy2 * 256 + idummy1
       byt2(1:1) = char(idummy1)
       return
       end subroutine num2bit2


* --------------------------------------------
* end of this file, thank you.
* --------------------------------------------
