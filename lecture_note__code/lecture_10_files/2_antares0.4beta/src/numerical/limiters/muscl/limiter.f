*--- MUSCL limiter
      real*8 function fslop(aa,bb,cc)
      real*8 aa,bb,cc
      fslop=(dsign(1.d0,(bb-aa))+dsign(1.d0,(cc-bb)))
     &  *dmin1(dmin1(0.5d0*dabs(cc-aa),
     &          2.d0*dabs(bb-aa) ),2.d0*dabs(cc-bb))   

      ! it seems has problem
      end