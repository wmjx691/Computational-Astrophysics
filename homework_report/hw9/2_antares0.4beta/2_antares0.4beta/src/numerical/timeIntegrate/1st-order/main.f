*===========================================================
*
* 1~3D Htdro and MHD codes      
*      
* Odie K.C. Pan, Hsiang-Hsu Wang, David C.C. Yen, and Chi Yuan
*      
*               Dec.   ,2007 1D and 2D
*               Mar.   ,2008 2D MHD
*               Mar.   ,2008 3D 
*               Mar. 15,2008 1~3D Hydro & MHD     
*               Apr. 13,2008 2D+3D Div clean 
*     

! 1st order in Time
!
!      
      program ANTARES
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'color.h'
      include 'io.h'
      include 'divFree.h'
      real*8 c1(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c2(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c3(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c4(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c5(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c6(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c7(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 c8(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8 uini(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf,1:8)
      real*8 cputime,elapsed
      real*8 clock(3)
      clock = 0
      t     = 0.d0
      tnext = 0.d0
      istep = 0
      nstop = 1000000000
      nprnt = 0
      dK    = (5.d0-3.d0*gam)/(gam-1.d0)

      nx = iend - ibeg + 1
      dx = (xmax-xmin)/(iend-ibeg+1)
      ny = jend - jbeg + 1
      dy = (ymax-ymin)/(jend-jbeg+1)
      nz = kend - kbeg + 1
      dz = (zmax-zmin)/(kend-kbeg+1)
!---  setupt the environment ------------------------------------------
      write(*,*) "====================================================="
      write(*,*) "  Problems Properity                                 "
      write(*,*) "====================================================="
      call checkconsistent() ! check consistent
      if (D1D)then
      write(*,*) "Dimension : 1D "
      elseif (D2D)then
      write(*,*) "Dimension : 2D "
      elseif (D3D)then
      write(*,*) "Dimension : 3D "
      else
      write(*,*) "error: please choose a dimension"
      stop      
      endif
      if (MHD)then
      write(*,*) "Gas       : MHD"
      else
      write(*,*) "Gas       : Hydro"
      endif
      if (ISOTHERMAL)then
              if(MHD) then 
                write(*,*)'error: currently, NO isothermal MHD !';stop
              endif   
      write(*,*) "E.O.S.    : Isothermal"
      else
      write(*,*) "E.O.S.    : Adiabatic"
      endif
      if (cartesian) then
      write(*,*) "geometry  : Cartesian"
      elseif(spherical)then
      write(*,*) "geometry  : Spherical"
      else
      write(*,*) "geometry  : cylindrical"
      endif
      write(*,*) "-----------------------------------------------------"
      write(*,*) "xmin=",xmin,"xmax=",xmax
      write(*,*) "ymin=",ymin,"ymax=",ymax
      write(*,*) "zmin=",zmin,"zmax=",zmax
      write(*,*) "dx=",dx,"dy=",dy,"dz=",dz
      write(*,*) "resolution:",nx,' X ',ny,' X ',nz
      if (ISOTHERMAL) then
      write(*,*) "isothermal sound speed =",iso_snd
      else
      write(*,*) "gamma=",gam
      endif
      write(*,*) "Simulation time=",tf
      write(*,*) "-----------------------------------------------------"
      write(*,*) "Setup complete"
      write(*,*) "-----------------------------------------------------"
      call grid


!---  SETUP initial condition for specific problem      
      call initial(den,px,py,pz,ene,bx,by,bz)

!---  STORE INITIAL CONDITION     
      do k=kbeg-kbuf,kend+kbuf
       do j=jbeg-jbuf,jend+jbuf
        do i=ibeg-ibuf,iend+ibuf
          uini(i,j,k,1) = den(i,j,k)
          uini(i,j,k,2) =  px(i,j,k)
          uini(i,j,k,3) =  py(i,j,k)
          uini(i,j,k,4) =  pz(i,j,k)
          if (.not. ISOTHERMAL) then
          uini(i,j,k,5) = ene(i,j,k)
          endif
          if(MHD)then
          uini(i,j,k,6) =  bx(i,j,k)
          uini(i,j,k,7) =  by(i,j,k)
          uini(i,j,k,8) =  bz(i,j,k)
          endif
        enddo
       enddo
      enddo
!--------------------

      call ddt(den,px,py,pz,ene,bx,by,bz,dx,dy,dt,dK)
      write(*,*)'initial dt =',dt 
      call io_initial(den,px,py,pz,ene,bx,by,bz,t,dt,color2,nprnt,UM)
!===================================================================
!--- START the MAIN LOOP

      do 1000 while (t.lt.tf.and.istep.le. nstop)
        call ddt(den,px,py,pz,ene,bx,by,bz,dx,dy,dt,dK)
        if (t+dt .gt. tf)then
             dt = tf-t
        endif
     
        
!=================================================================
*---  RK 1 --------------------------------------------------------
!=================================================================
!---  Boundary Condition        
        call boundary(den ,px ,py ,pz ,ene ,bx ,by ,bz,
     &                den ,px ,py ,pz ,ene ,bx ,by ,bz,uini )
!=================================================================
!                                REA 
!=================================================================
100     call godnov3d(t,dt,dx,1,
     &                 den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                 den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                 den2,px2,py2,pz2,ene2,bx2,by2,bz2,)
       if (D2D .or. D3D)then
        call godnov3d(t,dt,dy,2,
     &                 den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                 den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                 den2,px2,py2,pz2,ene2,bx2,by2,bz2)
       endif
       if (D3D)then
        call godnov3d(t,dt,dz,3,
     &                 den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                 den2,px2,py2,pz2,ene2,bx2,by2,bz2,
     &                 den2,px2,py2,pz2,ene2,bx2,by2,bz2)
       endif
!---  SOURCE TERM
        call source(den ,px ,py ,pz ,ene ,bx ,by ,bz,
     &              den2,px2,py2,pz2,ene2,bx2,by2,bz2,t,dt)
!---  DO DIVERGENCE CLEAN
      if (MHD.neqv.D1D .and. MHD) then 
      call divclean(dt,dx,dy,dz,
     &              den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &              den2,px2,py2,pz2,ene2,bx2,by2,bz2)
      endif



!---   checking values -------------------------------

        do k=kbeg,kend
        do j=jbeg,jend
        do i=ibeg,iend

         if (den2(i,j,k).le. 0.d0) then
           dt = 0.5d0*dt
           if(dt .le. 1.d-10) then
             write(*,*)'dt too small',dt
             stop
           endif
           write(*,*)'Oooops dt -> 0.5dt'
           goto 100
         endif
        enddo
        enddo
        enddo

        
        do k=kbeg,kend
        do j=jbeg,jend
        do i=ibeg,iend
          den(i,j,k)=den2(i,j,k)
           px(i,j,k)= px2(i,j,k)
           py(i,j,k)= py2(i,j,k)
           pz(i,j,k)= pz2(i,j,k)
           if (.not. ISOTHERMAL)then
            ene(i,j,k)=ene2(i,j,k)
           endif
           if(MHD)then
            bx(i,j,k)=bx2(i,j,k)
            by(i,j,k)=by2(i,j,k)
            bz(i,j,k)=bz2(i,j,k)
           endif
        enddo
        enddo
        enddo

      call io_loop(UM,den,px,py,pz,ene,bx,by,bz,t,dt,color2,nprnt)

      write(*,*)'t = ',t,'dt = ',dt,'istep = ',istep
      istep =istep + 1
      t =t +dt
1000  continue
!---  END OF MAIN LOOP 
!===================================================================

      write(*,*) "-----------------------------------------------------"
      write(*,*)'Calculation Finish !'
      write(*,*) "-----------------------------------------------------"
      call CPU_TIME(cputime)
        write(*,*)'CPU time',cputime,' sec'
        write(*,*)'Zone-cycles/sec',
     &  istep*(iend+4)*(jend+4)*(kend+4)/cputime
      write(*,*) "-----------------------------------------------------"
      stop 
      end program
