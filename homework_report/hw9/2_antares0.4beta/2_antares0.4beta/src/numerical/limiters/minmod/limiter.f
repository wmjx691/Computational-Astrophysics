*---  Minmod limiter      
      real*8 function fslop(aa,bb,cc)
      real*8 aa,bb,cc
      fslop=(dmax1(dsign(1.d0,(bb-aa)*(cc-bb)),0.d0)*
     &  dsign(1.d0,(cc-aa))*dmin1(dabs(0.5d0*(cc-aa)),
     &  dmin1(2.d0*dabs((bb-aa)),2.d0*dabs((cc-bb)))))
      end
