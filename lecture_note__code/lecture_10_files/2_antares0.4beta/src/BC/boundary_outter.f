      subroutine bnd_outter_outflow(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                             den2,px2,py2,pz2,ene2,bx2,by2,bz2,nd)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      integer nbc,ii,jj,kk,njbeg,njend,nkbeg,nkend
      integer n1,n2,n3,n4,n5,n6
      integer m1,m2,m3,m4,m5,m6
      integer a1,a2,a3,b1,b2,b3

      if (nd.eq.1)then
              nkbeg=kbeg
              nkend=kend
              njbeg=jbeg
              njend=jend
      elseif(nd.eq.2)then   
              nkbeg=kbeg
              nkend=kend
              njbeg=ibeg
              njend=iend
      elseif(nd.eq.3)then   
              nkbeg=jbeg
              nkend=jend
              njbeg=ibeg
              njend=iend
      endif

      do kk=nkbeg,nkend
      do jj=njbeg,njend
       if (nd .eq. 1)then
          n1=iend+ibuf
          n2=jj
          n3=kk
          n4=iend-1
          n5=jj
          n6=kk
          m1=iend+1
          m2=jj
          m3=kk
          m4=iend
          m5=jj
          m6=kk
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.2)then     
          n1=jj
          n2=jend+jbuf
          n3=kk
          n4=jj
          n5=jend-1
          n6=kk
          m1=jj
          m2=jend+1
          m3=kk
          m4=jj
          m5=jend
          m6=kk
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.3)then     
          n1=jj
          n2=kk
          n3=kend+kbuf
          n4=jj
          n5=kk
          n6=kend-1
          m1=jj
          m2=kk
          m3=kend+1
          m4=jj
          m5=kk
          m6=kend
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       endif
           
      den2(n1,n2,n3)= den(n4,n5,n6)
      den2(m1,m2,m3)= den(m4,m5,m6)
       px2(n1,n2,n3)= a1*px(n4,n5,n6)
       px2(m1,m2,m3)= a1*px(m4,m5,m6)
       py2(n1,n2,n3)= a2*py(n4,n5,n6)
       py2(m1,m2,m3)= a2*py(m4,m5,m6)
       pz2(n1,n2,n3)= a3*pz(n4,n5,n6)
       pz2(m1,m2,m3)= a3*pz(m4,m5,m6)
      if (.not. ISOTHERMAL)then
      ene2(n1,n2,n3)= ene(n4,n5,n6)
      ene2(m1,m2,m3)= ene(m4,m5,m6)
      endif
      if (MHD) then
       bx2(n1,n2,n3)= b1*bx(n4,n5,n6)
       bx2(m1,m2,m3)= b1*bx(m4,m5,m6)
       by2(n1,n2,n3)= b2*by(n4,n5,n6)
       by2(m1,m2,m3)= b2*by(m4,m5,m6)
       bz2(n1,n2,n3)= b3*bz(n4,n5,n6)
       bz2(m1,m2,m3)= b3*bz(m4,m5,m6)
      endif

      enddo
      enddo
      end
      subroutine bnd_outter_reflect(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                             den2,px2,py2,pz2,ene2,bx2,by2,bz2,nd)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      integer nbc,ii,jj,kk,njbeg,njend,nkbeg,nkend
      integer n1,n2,n3,n4,n5,n6
      integer m1,m2,m3,m4,m5,m6
      integer a1,a2,a3,b1,b2,b3

      if (nd.eq.1)then
              nkbeg=kbeg
              nkend=kend
              njbeg=jbeg
              njend=jend
      elseif(nd.eq.2)then   
              nkbeg=kbeg
              nkend=kend
              njbeg=ibeg
              njend=iend
      elseif(nd.eq.3)then   
              nkbeg=jbeg
              nkend=jend
              njbeg=ibeg
              njend=iend
      endif

      do kk=nkbeg,nkend
      do jj=njbeg,njend
       if (nd .eq. 1)then
          n1=iend+ibuf
          n2=jj
          n3=kk
          n4=iend-1
          n5=jj
          n6=kk
          m1=iend+1
          m2=jj
          m3=kk
          m4=iend
          m5=jj
          m6=kk
          a1=-1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.2)then     
          n1=jj
          n2=jend+jbuf
          n3=kk
          n4=jj
          n5=jend-1
          n6=kk
          m1=jj
          m2=jend+1
          m3=kk
          m4=jj
          m5=jend
          m6=kk
          a1= 1.d0
          a2=-1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.3)then     
          n1=jj
          n2=kk
          n3=kend+kbuf
          n4=jj
          n5=kk
          n6=kend-1
          m1=jj
          m2=kk
          m3=kend+1
          m4=jj
          m5=kk
          m6=kend
          a1= 1.d0
          a2= 1.d0
          a3=-1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       endif
           
      den2(n1,n2,n3)= den(n4,n5,n6)
      den2(m1,m2,m3)= den(m4,m5,m6)
       px2(n1,n2,n3)= a1*px(n4,n5,n6)
       px2(m1,m2,m3)= a1*px(m4,m5,m6)
       py2(n1,n2,n3)= a2*py(n4,n5,n6)
       py2(m1,m2,m3)= a2*py(m4,m5,m6)
       pz2(n1,n2,n3)= a3*pz(n4,n5,n6)
       pz2(m1,m2,m3)= a3*pz(m4,m5,m6)
      if (.not. ISOTHERMAL)then
      ene2(n1,n2,n3)= ene(n4,n5,n6)
      ene2(m1,m2,m3)= ene(m4,m5,m6)
      endif
      if (MHD) then
       bx2(n1,n2,n3)= b1*bx(n4,n5,n6)
       bx2(m1,m2,m3)= b1*bx(m4,m5,m6)
       by2(n1,n2,n3)= b2*by(n4,n5,n6)
       by2(m1,m2,m3)= b2*by(m4,m5,m6)
       bz2(n1,n2,n3)= b3*bz(n4,n5,n6)
       bz2(m1,m2,m3)= b3*bz(m4,m5,m6)
      endif

      enddo
      enddo
      end
      subroutine bnd_outter_periodic(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                             den2,px2,py2,pz2,ene2,bx2,by2,bz2,nd)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      integer nbc,ii,jj,kk,njbeg,njend,nkbeg,nkend
      integer n1,n2,n3,n4,n5,n6
      integer m1,m2,m3,m4,m5,m6
      integer a1,a2,a3,b1,b2,b3

      if (nd.eq.1)then
              nkbeg=kbeg
              nkend=kend
              njbeg=jbeg
              njend=jend
      elseif(nd.eq.2)then   
              nkbeg=kbeg
              nkend=kend
              njbeg=ibeg
              njend=iend
      elseif(nd.eq.3)then   
              nkbeg=jbeg
              nkend=jend
              njbeg=ibeg
              njend=iend
      endif

      do kk=nkbeg,nkend
      do jj=njbeg,njend
       if (nd .eq. 1)then
          n1=iend+ibuf
          n2=jj
          n3=kk
          n4=ibeg+1
          n5=jj
          n6=kk
          m1=iend+1
          m2=jj
          m3=kk
          m4=ibeg
          m5=jj
          m6=kk
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.2)then     
          n1=jj
          n2=jend+jbuf
          n3=kk
          n4=jj
          n5=jbeg+1
          n6=kk
          m1=jj
          m2=jend+1
          m3=kk
          m4=jj
          m5=jbeg
          m6=kk
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       elseif(nd.eq.3)then     
          n1=jj
          n2=kk
          n3=kend+kbuf
          n4=jj
          n5=kk
          n6=kbeg+1
          m1=jj
          m2=kk
          m3=kend+1
          m4=jj
          m5=kk
          m6=kbeg
          a1= 1.d0
          a2= 1.d0
          a3= 1.d0
          b1= 1.d0
          b2= 1.d0
          b3= 1.d0
       endif
           
      den2(n1,n2,n3)= den(n4,n5,n6)
      den2(m1,m2,m3)= den(m4,m5,m6)
       px2(n1,n2,n3)= a1*px(n4,n5,n6)
       px2(m1,m2,m3)= a1*px(m4,m5,m6)
       py2(n1,n2,n3)= a2*py(n4,n5,n6)
       py2(m1,m2,m3)= a2*py(m4,m5,m6)
       pz2(n1,n2,n3)= a3*pz(n4,n5,n6)
       pz2(m1,m2,m3)= a3*pz(m4,m5,m6)
      if (.not. ISOTHERMAL)then
      ene2(n1,n2,n3)= ene(n4,n5,n6)
      ene2(m1,m2,m3)= ene(m4,m5,m6)
      endif
      if (MHD) then
       bx2(n1,n2,n3)= b1*bx(n4,n5,n6)
       bx2(m1,m2,m3)= b1*bx(m4,m5,m6)
       by2(n1,n2,n3)= b2*by(n4,n5,n6)
       by2(m1,m2,m3)= b2*by(m4,m5,m6)
       bz2(n1,n2,n3)= b3*bz(n4,n5,n6)
       bz2(m1,m2,m3)= b3*bz(m4,m5,m6)
      endif

      enddo
      enddo
      end
      subroutine bnd_outter_ini(den ,px ,py ,pz ,ene ,bx ,by ,bz ,
     &                      den2,px2,py2,pz2,ene2,bx2,by2,bz2,nd,uini)
      implicit none
      include 'fluid.h'
      include 'parameters'
      include 'conserve.h'
      integer nbc,ii,jj,kk,njbeg,njend,nkbeg,nkend
      integer n1,n2,n3,n4,n5,n6
      integer m1,m2,m3,m4,m5,m6
      integer a1,a2,a3,b1,b2,b3
      real*8  uini(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf,
     &             kbeg-kbuf:kend+kbuf,1:8)

      if (nd.eq.1)then
              nkbeg=kbeg
              nkend=kend
              njbeg=jbeg
              njend=jend
      elseif(nd.eq.2)then   
              nkbeg=kbeg
              nkend=kend
              njbeg=ibeg
              njend=iend
      elseif(nd.eq.3)then   
              nkbeg=jbeg
              nkend=jend
              njbeg=ibeg
              njend=iend
      endif

      do kk=nkbeg,nkend
      do jj=njbeg,njend
       if (nd .eq. 1)then
          n1=iend+ibuf
          n2=jj
          n3=kk
          n4=iend+ibuf
          n5=jj
          n6=kk
          m1=iend+1
          m2=jj
          m3=kk
          m4=iend+1
          m5=jj
          m6=kk
       elseif(nd.eq.2)then     
          n1=jj
          n2=jend+jbuf
          n3=kk
          n4=jj
          n5=jend+jbuf
          n6=kk
          m1=jj
          m2=jend+1
          m3=kk
          m4=jj
          m5=jend+1
          m6=kk
       elseif(nd.eq.3)then     
          n1=jj
          n2=kk
          n3=kend+kbuf
          n4=jj
          n5=kk
          n6=kend+kbuf
          m1=jj
          m2=kk
          m3=kend+1
          m4=jj
          m5=kk
          m6=kend+1
       endif
      den2(n1,n2,n3)= uini(n4,n5,n6,1)
      den2(m1,m2,m3)= uini(m4,m5,m6,1)
       px2(n1,n2,n3)= uini(n4,n5,n6,2)
       px2(m1,m2,m3)= uini(m4,m5,m6,2)
       py2(n1,n2,n3)= uini(n4,n5,n6,3)
       py2(m1,m2,m3)= uini(m4,m5,m6,3)
       pz2(n1,n2,n3)= uini(n4,n5,n6,4)
       pz2(m1,m2,m3)= uini(m4,m5,m6,4)
      if (.not. ISOTHERMAL)then
      ene2(n1,n2,n3)= uini(n4,n5,n6,5)
      ene2(m1,m2,m3)= uini(m4,m5,m6,5)
      endif
      if (MHD) then
       bx2(n1,n2,n3)= uini(n4,n5,n6,6)
       bx2(m1,m2,m3)= uini(m4,m5,m6,6)
       by2(n1,n2,n3)= uini(n4,n5,n6,7)
       by2(m1,m2,m3)= uini(m4,m5,m6,7)
       bz2(n1,n2,n3)= uini(n4,n5,n6,8)
       bz2(m1,m2,m3)= uini(m4,m5,m6,8)
      endif

      enddo
      enddo
      end
