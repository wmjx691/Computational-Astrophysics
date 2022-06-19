      subroutine io_initial(den,px,py,pz,ene,bx,by,bz,
     &                      t,dt,color2,nprnt,UM)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'io.h'
      include 'color.h'
      call colormap(color2)

      if (fio_bin_o .gt. fio_ppm) then
              write(*,*)'please set fio_bin less than fio_ppm'
      endif
      if (fio_bin_o.ge. 1.d0 ) then
        dprnt_bin = tf/(fio_bin_o+1.d-10)
      else
        dprnt_bin = 0.d0
      endif      
      if (fio_ppm.ge. 1.d0 ) then
        dprnt_ppm = tf/(fio_ppm+1.d-10)
      else
        dprnt_ppm = 0.d0
      endif      
      nprnt = 0
      dnext_bin = dprnt_bin
      dnext_ppm = dprnt_ppm
*---  output or reset the initial data
      if (fio_bin_i .eq. 0) then
              if(.not.D1D)then
              call io_ppm(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0,color2)
              endif
              call io_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0)
      else
              nprnt = fio_bin_i
              call io_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,1)
              dnext_bin = t+dprnt_bin
              dnext_ppm = t+dprnt_ppm
      endif
      end

      subroutine io_loop(UM,den,px,py,pz,ene,bx,by,bz,t,dt,color2,nprnt)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      include 'io.h'
      include 'color.h'
      if (fio_ppm.eq. 0.d0 )then
              nprnt = 1 + nprnt
              if(.not.D1D)then
              call io_ppm(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0,color2)
              endif
              if(fio_bin_o .eq. 0.d0) then
                call io_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0)  
              elseif(t.gt. dnext_bin.or. t.ge.tf) then 
                dnext_bin = dnext_bin + dprnt_bin
                call io_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0)  
              endif  
      elseif(t.gt. dnext_ppm .or. t.ge.tf) then
              nprnt = 1 + nprnt
              dnext_ppm = dnext_ppm + dprnt_ppm      
              if(.not.D1D)then
              call io_ppm(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0,color2)
              endif
              if (t.gt.dnext_bin .or. t.ge.tf) then            
              dnext_bin = dnext_bin + dprnt_bin
              call io_bin(nprnt,t,UM,den,px,py,pz,ene,bx,by,bz,0)  
              endif
      endif
      end
