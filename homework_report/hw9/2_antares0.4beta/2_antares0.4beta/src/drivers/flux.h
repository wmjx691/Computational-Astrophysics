      integer m,cx,cy,cz,f1,f2,f3,f4,f5,f6,f7,f8
      real*8  g1,g2,g3,g4,g5,g6,g7,g8
      real*8  aa,bb,cc
      real*8  coef(8),UL(8),UR(8),slopeL(8),slopeR(8)
      real*8  F(8),VL(8),VR(8)	
      real*8  dtdx,snd
      real*8  rhoL,rhoM,rhoR
      real*8  vxL,vxM,vxR,pxL,pxM,pxR    
      real*8  vyL,vyM,vyR,pyL,pyM,pyR    
      real*8  vzL,vzM,vzR,pzL,pzM,pzR 
      real*8  enL,enM,enR 
      real*8  bxL,bxM,bxR    
      real*8  byL,byM,byR    
      real*8  bzL,bzM,bzR    
      real*8  S_den(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf  
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_vx(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_vy(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf  
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_vz(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_ene(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_bx(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_by(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  S_bz(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  SL(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf,1:8)
      real*8  F_den(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_px(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_py(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_pz(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_ene(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_bx(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_by(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)
      real*8  F_bz(ibeg-ibuf:iend+ibuf,jbeg-jbuf:jend+jbuf
     &          ,kbeg-kbuf:kend+kbuf)