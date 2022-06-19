      integer IOC
      real*8 dnext_ppm,dnext_bin,tnext
      real*8 dprnt_ppm,dprnt_bin
      real*8 UM(ibeg:iend+2*ibuf,jbeg:jend+2*jbuf
     &         ,kbeg:kend+2*kbuf,1:8)
      common /ppm/ dnext_ppm,dprnt_ppm
      common /bin/ dnext_bin,dprnt_bin
