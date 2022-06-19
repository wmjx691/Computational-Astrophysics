*--- Van Leer limiter
      real*8 function fslop(aa,bb,cc)
      real*8 aa,bb,cc
      fslop=((dsign(1.d0,(bb-aa))+dsign(1.d0,(cc-bb)))
     &  *(dabs((bb-aa))*dabs((cc-bb)))/
     &   (dabs((bb-aa))+dabs((cc-bb))+1.d-8))
      end
