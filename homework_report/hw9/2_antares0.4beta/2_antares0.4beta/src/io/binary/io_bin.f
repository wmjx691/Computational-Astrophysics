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

       subroutine IO_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,IOC)
       implicit none
       include 'fluid.h'
       include 'parameters'
       include 'conserve.h'
       include 'io.h'
       include 'color.h'
       integer   ihpixf,jvpixf,khpixf

       ihpixf = iend+2*ibuf
       jvpixf = jend+2*jbuf       
       khpixf = kend+2*kbuf
       if (D1D) jvpixf=1
       if (D2D) khpixf=1

       call  binout(den,px,py,pz,ene,bx,by,bz,
     &              UM,t,nprnt,ihpixf,jvpixf,khpixf,IOC)
       end

       subroutine binout(den,px,py,pz,ene,bx,by,bz,
     &                   UM,t,nframe,ihpixf,jvpixf,khpixf,IOC)
       implicit none
       include 'fluid.h'
       include 'parameters'
       include 'conserve.h'
       include 'io.h'
       integer ihpixf, jvpixf, khpixf,n
       real*8 Uhydro(ihpixf,jvpixf,khpixf,5)
       real*8   Umhd(ihpixf,jvpixf,khpixf,8)
       integer nframe
       character*12 fnameout
       integer idummy, icnt,iomethod
       real*8  rdummy
       character*14 frmtstr
       character*54 headmsw
       character*4  byt4
       character*2  byt2
       real*8  cc
       integer*8 sizeofreal
       parameter (sizeofreal=8)

       if (ioc .EQ. 0) then
!--- 1D
         if (D1D)then
         write(fnameout,'(''data'',i4.4,''.tab'')') nframe ! name of bin file
         open(unit=3,file=fnameout,status='unknown')
         write(*,*) 'Now writing tab      file : ', fnameout
         j=1
         k=1
          do i=ibeg,iend
            write(3,*)i,x(i),den(i,j,k),px(i,j,k)/den(i,j,k),ene(i,j,k)
          enddo
         close(3)
!--- 2D
         elseif(D2D)then

        k=1
        do j=1,jvpixf
        do i=1,ihpixf
                if (.not. MHD) then
                Uhydro(i,j,k,1) = den(i-2,j-2,k)
                Uhydro(i,j,k,2) =  px(i-2,j-2,k)
                Uhydro(i,j,k,3) =  py(i-2,j-2,k)
                Uhydro(i,j,k,4) =  pz(i-2,j-2,k)
                Uhydro(i,j,k,5) = ene(i-2,j-2,k)
                else
                Umhd(i,j,k,1) = den(i-2,j-2,k)
                Umhd(i,j,k,2) =  px(i-2,j-2,k)
                Umhd(i,j,k,3) =  py(i-2,j-2,k)
                Umhd(i,j,k,4) =  pz(i-2,j-2,k)
                Umhd(i,j,k,5) = ene(i-2,j-2,k)
                Umhd(i,j,k,6) =  bx(i-2,j-2,k)
                Umhd(i,j,k,7) =  by(i-2,j-2,k)
                Umhd(i,j,k,8) =  bz(i-2,j-2,k)
                endif
        enddo
        enddo

         write(fnameout,'(''dat'',i4.4,''.bin'')') nframe ! name of bin file
         !open(unit=60,file=fnameout,status='unknown',form='unformatted')
         open(unit=60,file=fnameout,form='unformatted')
         write(*,*) 'Now writing tab      file : ', fnameout

         if (.not. MHD) then

        ! write(60) ihpixf,jvpixf,khpixf,5
         ! do n=1,5
         ! do k=1,khpixf
         ! do j=1,jvpixf
         !    write(60) (Uhydro(i,j,k,n),i=1,ihpixf)
         ! enddo
         ! enddo
         ! enddo 
           write(60) Uhydro      
         close(60)

         else

         !write(60) ihpixf,jvpixf,khpixf,8
         ! do n=1,8
         ! do k=1,khpixf
         ! do j=1,jvpixf
         !   write(60) (Umhd(i,j,k,n),i=1,ihpixf)
         ! enddo
         ! enddo
         ! enddo 
          write(60) Umhd      
         close(60)
         endif
!--- 3D
         else        
* bin
        !change den,px,... to U


        do k=1,khpixf
        do j=1,jvpixf
        do i=1,ihpixf
                if (.not. MHD) then
                Uhydro(i,j,k,1) = den(i-2,j-2,k-2)
                Uhydro(i,j,k,2) =  px(i-2,j-2,k-2)
                Uhydro(i,j,k,3) =  py(i-2,j-2,k-2)
                Uhydro(i,j,k,4) =  pz(i-2,j-2,k-2)
                Uhydro(i,j,k,5) = ene(i-2,j-2,k-2)
                else
                Umhd(i,j,k,1) = den(i-2,j-2,k-2)
                Umhd(i,j,k,2) =  px(i-2,j-2,k-2)
                Umhd(i,j,k,3) =  py(i-2,j-2,k-2)
                Umhd(i,j,k,4) =  pz(i-2,j-2,k-2)
                Umhd(i,j,k,5) = ene(i-2,j-2,k-2)
                Umhd(i,j,k,6) =  bx(i-2,j-2,k-2)
                Umhd(i,j,k,7) =  by(i-2,j-2,k-2)
                Umhd(i,j,k,8) =  bz(i-2,j-2,k-2)
                endif
        enddo
        enddo
        enddo

         write(fnameout,'(''dat'',i4.4,''.bin'')') nframe ! name of bin file
         !open(unit=60,file=fnameout,status='unknown',form='unformatted')
         open(unit=1,file=fnameout,form='unformatted')
         write(*,*) 'Now writing tab      file : ', fnameout

         if (.not. MHD) then

          !write(60) ihpixf,jvpixf,khpixf,5
          !do n=1,5
          !do k=1,khpixf
          !do j=1,jvpixf
          !   write(60) (Uhydro(i,j,k,n),i=1,ihpixf)
          !enddo
          !enddo
          !enddo       

          write(1) Uhydro
         close(1)

         else

         !write(60) ihpixf,jvpixf,khpixf,8
         ! do n=1,8
         ! do k=1,khpixf
         ! do j=1,jvpixf
         !   write(60) (Umhd(i,j,k,n),i=1,ihpixf)
         ! enddo
         ! enddo
         ! enddo 
         write(1) Umhd      
         close(1)
         endif

         endif ! Endif D2D,D3D
*=============================================
* for input
*
*
*--------------------------------------------       
       elseif (ioc .EQ. 1) then

         write(*,*)'Error: data input no support yet.'
         stop 
       else 
         write(*,*)'Error: no such ioc'      
         stop 
       endif

       return
       end 