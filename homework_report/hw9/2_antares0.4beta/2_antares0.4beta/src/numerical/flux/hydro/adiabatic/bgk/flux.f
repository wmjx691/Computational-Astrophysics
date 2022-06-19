c    "A Gas-Kinetic BGK Scheme for the Navier-Stokes Equations
c    and Its connection with Artificial Dissipation and Godunov
c    Method", Kun Xu, Journal of Computational Physics, No 171,
c    289-335 (2001)

ccccccc    invM.f is called in this subroutine  

ccccccc    Hsiang-Hsu Wang Jan. 2007
      !subroutine flux(WL,WR,slopeL,slopeR,dx,dt,K,FF)
      subroutine flux(WL,WR,slopeL,slopeR,FF,dt,sndL,sndR)
      include 'fluid.h'
      include 'parameters'
      dimension WL(5), WR(5), slopeL(5), slopeR(5), F(5),FF(5)
      dimension gamma1_mul(4), gamma2_mul(4), gamma3_mul(4)
      dimension gamma4_mul(4), gamma5_mul(4)
      dimension beta1_mul(4), beta2_mul(4)
      dimension beta3_mul(4), beta4_mul(4)
      dimension temp1_mul(4), temp2_mul(4)
      real*8 cK
      real*8 WL, WR, slopeL, slopeR, F,FF
      real*8 gamma1_mul, gamma2_mul, gamma3_mul
      real*8 gamma4_mul, gamma5_mul
      real*8 beta1_mul, beta2_mul, beta3_mul, beta4_mul
      real*8 temp1_mul, temp2_mul
      real*8 pi, rhoL, XmomL, YmomL, engL, UL, VL
      real*8     rhoR, XmomR, YmomR, engR, UR, VR
      real*8 lamdaL, lamdaR
      real*8 a1L, a2L, a3L, a4L, a1R, a2R, a3R, a4R
      real*8 xi2L, xi4L, xi2R, xi4R
      real*8 U0L, U1L, U2L, U3L, U4L, U5L, U6L
      real*8 V0L, V1L, V2L, V3L, V4L, V5L, V6L
      real*8 U0R, U1R, U2R, U3R, U4R, U5R, U6R
      real*8 V0R, V1R, V2R, V3R, V4R, V5R, V6R
      real*8 ALb1, ALb2, ALb3, ALb4, ARb1, ARb2, ARb3, ARb4
      real*8 AL1, AL2, AL3, AL4, AR1, AR2, AR3, AR4
      real*8 U0GL, U1GL, U2GL, U3GL, U4GL, U5GL, U6GL
      real*8 U0LR, U1LR, U2LR, U3LR, U4LR, U5LR, U6LR
      real*8 rho0, Xmom0, Ymom0, eng0, lamda0, xi2, xi4, U0, V0
      real*8 abarL_b1, abarL_b2, abarL_b3, abarL_b4
      real*8 abarR_b1, abarR_b2, abarR_b3, abarR_b4
      real*8 a1barL, a2barL, a3barL, a4barL
      real*8 a1barR, a2barR, a3barR, a4barR
      real*8 appL, appR, tau, eta, eta1
      real*8 gamma1, gamma2, gamma3, gamma4, gamma5
      real*8 U0_0, U1_0, U2_0, U3_0, U4_0, U5_0, U6_0
      real*8 V0_0, V1_0, V2_0, V3_0, V4_0, V5_0, V6_0
      real*8 U0G_0,U1G_0,U2G_0,U3G_0,U4G_0,U5G_0,U6G_0
      real*8 U0L_0,U1L_0,U2L_0,U3L_0,U4L_0,U5L_0,U6L_0
      real*8 Ab1,Ab2,Ab3,Ab4,A1,A2,A3,A4 
      real*8 beta1,beta2, beta3, beta4, beta5

      write(*,*)'warnning: BGK solver has not modified for this'
      write(*,*)'          version of antares codes.'
      write(*,*)'          try different flux solver.'
      stop

      if (D3D) then
        write(*,*)'Error: BGK solver only work for 1D or 2D'
        stop
      endif 

      dx = x(2)-x(1)
      cK  = (4.d0-2.d0*gam)/(gam-1.d0)
      pi = 4.d0*datan(1.d0)

      rhoL  = WL(1) + 0.5d0*slopeL(1)*dx
      XmomL = (WL(2) + 0.5d0*slopeL(2)*dx)*rhoL
      YmomL = (WL(3) + 0.5d0*slopeL(3)*dx)*rhoL
      engL  = WL(5) + 0.5d0*slopeL(5)*dx

      rhoR  = WR(1) - 0.5d0*slopeR(1)*dx
      XmomR = (WR(2) - 0.5d0*slopeR(2)*dx)*rhoR
      YmomR = (WR(3) - 0.5d0*slopeR(3)*dx)*rhoR
      engR  = WR(5) - 0.5d0*slopeR(5)*dx


      UL = XmomL/rhoL
      VL = YmomL/rhoL
      
      UR = XmomR/rhoR
      VR = YmomR/rhoR


      lamdaL = 0.25d0*(cK+2.d0)*rhoL/(engL-0.5d0*(XmomL*XmomL 
     + +YmomL*YmomL)/rhoL)    
      lamdaR = 0.25d0*(cK+2.d0)*rhoR/(engR-0.5d0*(XmomR*XmomR 
     + +YmomR*YmomR)/rhoR)


      if(lamdaL .le. 0.d0 .or. lamdaR .le. 0.d0) then
      rhoL  = WL(1)
      XmomL = WL(2)*WL(1)
      YmomL = WL(3)*WL(1)
      engL  = WL(5)
      slopeL(1) = 0.d0
      slopeL(2) = 0.d0
      slopeL(3) = 0.d0
      slopeL(5) = 0.d0

      lamdaL = 0.25d0*(cK+2.d0)*rhoL/(engL-0.5d0*(XmomL*XmomL 
     + +YmomL*YmomL)/rhoL)
      UL = XmomL/rhoL
      VL = YmomL/rhoL

      rhoR  = WR(1)
      XmomR = WR(2)*WR(1)
      YmomR = WR(3)*WR(1)
      engR  = WR(5)
      slopeR(1) = 0.d0
      slopeR(2) = 0.d0
      slopeR(3) = 0.d0
      slopeR(5) = 0.d0

      lamdaR = 0.25d0*(cK+2.d0)*rhoR/(engR-0.5d0*(XmomR*XmomR 
     + + YmomR*YmomR)/rhoR) 
     
      UR = XmomR/rhoR
      VR = YmomR/rhoR 
      write(*,*) 'oops'
      endif
ccccccc (a1L,a2L,a3L,a4L)*rhoL  && (a1R,a2R,a3R,a4R)*rhoR
c      write(*,*) 'hello', lamdaL, lamdaR
      call invM(UL,VL,lamdaL,cK,slopeL(1),slopeL(2),slopeL(3),slopeL(5),
     +a1L,a2L,a3L,a4L)
      call invM(UR,VR,lamdaR,cK,slopeR(1),slopeR(2),slopeR(3),slopeR(5),
     +a1R,a2R,a3R,a4R)
cccccccccccccccccccccccccccccccccccccccccccccccc
      xi2L = 0.5d0*cK/lamdaL
      xi4L = 0.25d0*cK*(cK+2.d0)/(lamdaL*lamdaL)

      xi2R = 0.5d0*cK/lamdaR
      xi4R = 0.25d0*cK*(cK+2.d0)/(lamdaR*lamdaR) 

      U0L = 1.d0
      U1L = UL
      U2L = UL*U1L + 0.5d0/lamdaL*U0L
      U3L = UL*U2L + 1.d0/lamdaL*U1L
      U4L = UL*U3L + 1.5d0/lamdaL*U2L
      U5L = UL*U4L + 2.d0/lamdaL*U3L
      U6L = UL*U5L + 2.5d0/lamdaL*U4L
      
      V0L = 1.d0
      V1L = VL
      V2L = VL*V1L + 0.5d0/lamdaL*V0L
      V3L = VL*V2L + 1.d0/lamdaL*V1L
      V4L = VL*V3L + 1.5d0/lamdaL*V2L
      V5L = VL*V4L + 2.d0/lamdaL*V3L
      V6L = VL*V5L + 2.5d0/lamdaL*V4L      

      U0R = 1.d0
      U1R = UR
      U2R = UR*U1R + 0.5d0/lamdaR*U0R
      U3R = UR*U2R + 1.d0/lamdaR*U1R
      U4R = UR*U3R + 1.5d0/lamdaR*U2R
      U5R = UR*U4R + 2.d0/lamdaR*U3R
      U6R = UR*U5R + 2.5d0/lamdaR*U4R

      V0R = 1.d0
      V1R = VR
      V2R = VR*V1R + 0.5d0/lamdaR*V0R
      V3R = VR*V2R + 1.d0/lamdaR*V1R
      V4R = VR*V3R + 1.5d0/lamdaR*V2R
      V5R = VR*V4R + 2.d0/lamdaR*V3R
      V6R = VR*V5R + 2.5d0/lamdaR*V4R

      ALb1 = a1L*U1L + a2L*U2L + a3L*U1L*V1L + 
     +       0.5d0*a4L*(U3L+U1L*V2L+xi2L*U1L)
      ALb2 = a1L*U2L + a2L*U3L + a3L*U2L*V1L + 
     +       0.5d0*a4L*(U4L+U2L*V2L+xi2L*U2L)
      ALb3 = a1L*U1L*V1L + a2L*U2L*V1L + a3L*U1L*V2L +
     +       0.5d0*a4L*(U3L*V1L+U1L*V3L+xi2L*U1L*V1L)
      ALb4 = a1L*0.5d0*(U3L+U1L*V2L+U1L*xi2L) + 
     +       a2L*0.5d0*(U4L+U2L*V2L+U2L*xi2L) +
     +       a3L*0.5d0*(U3L*V1L+U1L*V3L+U1L*V1L*xi2L) +
     +       a4L*0.25d0*(U5L+U1L*V4L+U1L*xi4L+2.d0*U3L*V2L+
     +       2.d0*U3L*xi2L+2.d0*U1L*V2L*xi2L)
      
      ARb1 = a1R*U1R + a2R*U2R + a3R*U1R*V1R +
     +       0.5d0*a4R*(U3R+U1R*V2R+xi2R*U1R)
      ARb2 = a1R*U2R + a2R*U3R + a3R*U2R*V1R +
     +       0.5d0*a4R*(U4R+U2R*V2R+xi2R*U2R)
      ARb3 = a1R*U1R*V1R + a2R*U2R*V1R + a3R*U1R*V2R +
     +       0.5d0*a4R*(U3R*V1R+U1R*V3R+xi2R*U1R*V1R)
      ARb4 = a1R*0.5d0*(U3R+U1R*V2R+U1R*xi2R) +
     +       a2R*0.5d0*(U4R+U2R*V2R+U2R*xi2R) +
     +       a3R*0.5d0*(U3R*V1R+U1R*V3R+U1R*V1R*xi2R) +
     +       a4R*0.25d0*(U5R+U1R*V4R+U1R*xi4R+2.d0*U3R*V2R
     +       +2.d0*U3R*xi2R+2.d0*U1R*V2R*xi2R)


cccccc  (AL1,AL2,AL3,AL4)*rhoL && (AR1,AR2,AR3,AR4)*rhoR
      call invM(UL,VL,lamdaL,cK,-ALb1,-ALb2,-ALb3,-ALb4,AL1,AL2,AL3,AL4)
      call invM(UR,VR,lamdaR,cK,-ARb1,-ARb2,-ARb3,-ARb4,AR1,AR2,AR3,AR4)
ccccccccccccccccccccccc
c      write(*,*) WL(1), WL(2), WL(3), WL(4), lamdaL

      U0GL = 0.5d0*derfc(-dsqrt(lamdaL)*UL)
      U1GL = UL*U0GL + 0.5d0*dexp(-lamdaL*UL*UL)/dsqrt(pi*lamdaL)
      U2GL = UL*U1GL + 0.5d0*U0GL/lamdaL
      U3GL = UL*U2GL + 1.d0*U1GL/lamdaL
      U4GL = UL*U3GL + 1.5d0*U2GL/lamdaL
      U5GL = UL*U4GL + 2.d0*U3GL/lamdaL
      U6GL = UL*U5GL + 2.5d0*U4GL/lamdaL

      U0LR = 0.5d0*derfc(dsqrt(lamdaR)*UR)
      U1LR = UR*U0LR - 0.5d0*dexp(-lamdaR*UR*UR)/dsqrt(pi*lamdaR)
      U2LR = UR*U1LR + 0.5d0*U0LR/lamdaR
      U3LR = UR*U2LR + 1.d0*U1LR/lamdaR
      U4LR = UR*U3LR + 1.5d0*U2LR/lamdaR
      U5LR = UR*U4LR + 2.d0*U3LR/lamdaR
      U6LR = UR*U5LR + 2.5d0*U4LR/lamdaR
ccccccccccccc  rho0, Xmom0, Ymom0, lamda0, U0, V0
      rho0 = rhoL*U0GL + rhoR*U0LR
      Xmom0 = rhoL*U1GL + rhoR*U1LR
      Ymom0 = rhoL*U0GL*V1L + rhoR*U0LR*V1R
      eng0 = 0.5d0*(rhoL*(U2GL+U0GL*V2L+U0GL*xi2L)
     +           +rhoR*(U2LR+U0LR*V2R+U0LR*xi2R))

      
      lamda0 = 0.25d0*(cK+2.d0)*rho0/(eng0
     +           -0.5d0*(Xmom0*Xmom0+Ymom0*Ymom0)/rho0)  
      U0 = Xmom0/rho0
      V0 = Ymom0/rho0
      xi2 = 0.5d0*cK/lamda0
      xi4 = 0.25d0*cK*(cK+2.d0)/(lamda0*lamda0)

c      write(*,*) K
cccccccccccccc (a1barL/R, a2barL/R, a3barL/R)*rho0  
      abarL_b1 = 2.d0*(rho0-WL(1))/dx
      abarL_b2 = 2.d0*(Xmom0-WL(2))/dx
      abarL_b3 = 2.d0*(Ymom0-WL(3))/dx
      abarL_b4 = 2.d0*(eng0-WL(5))/dx

      abarR_b1 = 2.d0*(WR(1)-rho0)/dx
      abarR_b2 = 2.d0*(WR(2)-Xmom0)/dx
      abarR_b3 = 2.d0*(WR(3)-Ymom0)/dx
      abarR_b4 = 2.d0*(WR(5)-eng0)/dx
    
      call invM(U0,V0,lamda0,cK,abarL_b1,abarL_b2,abarL_b3,abarL_b4,
     +          a1barL,a2barL,a3barL,a4barL)
      call invM(U0,V0,lamda0,cK,abarR_b1,abarR_b2,abarR_b3,abarR_b4,
     +          a1barR,a2barR,a3barR,a4barR)
cccccccccccccccccccccc  (A1,A2,A3)*rho0  
      appL = rhoL/(2.d0*lamdaL)
      appR = rhoR/(2.d0*lamdaR)
       
      tau = 0.05d0*dt+dt*dabs(appL-appR)/(appL+appR)
      eta = dexp(-dt/tau)
      eta1 = 1.d0-eta

      gamma0 = dt-tau*eta1
      gamma1 = -eta1/gamma0
      gamma2 = (-dt + 2.d0*tau*eta1-dt*eta)/gamma0
      gamma3 = eta1/gamma0
      gamma4 = (dt*eta - tau*eta1)/gamma0
      gamma5 = tau*eta1/gamma0

      U0_0 = 1.d0
      U1_0 = U0
      U2_0 = U0*U1_0 + 0.5d0/lamda0*U0_0
      U3_0 = U0*U2_0 + 1.d0/lamda0*U1_0
      U4_0 = U0*U3_0 + 1.5d0/lamda0*U2_0
      U5_0 = U0*U4_0 + 2.d0/lamda0*U3_0
      U6_0 = U0*U5_0 + 2.5d0/lamda0*U4_0
 
      V0_0 = 1.d0
      V1_0 = V0
      V2_0 = V0*V1_0 + 0.5d0/lamda0*V0_0
      V3_0 = V0*V2_0 + 1.d0/lamda0*V1_0
      V4_0 = V0*V3_0 + 1.5d0/lamda0*V2_0
      V5_0 = V0*V4_0 + 2.d0/lamda0*V3_0
      V6_0 = V0*V5_0 + 2.5d0/lamda0*V4_0

      U0G_0 = 0.5d0*derfc(-dsqrt(lamda0)*U0)
      U1G_0 = U0*U0G_0 + 0.5d0*dexp(-lamda0*U0*U0)/dsqrt(pi*lamda0)
      U2G_0 = U0*U1G_0 + 0.5d0*U0G_0/lamda0
      U3G_0 = U0*U2G_0 + 1.d0*U1G_0/lamda0
      U4G_0 = U0*U3G_0 + 1.5d0*U2G_0/lamda0
      U5G_0 = U0*U4G_0 + 2.d0*U3G_0/lamda0
      U6G_0 = U0*U5G_0 + 2.5d0*U4G_0/lamda0

      U0L_0 = 0.5d0*derfc(dsqrt(lamda0)*U0)
      U1L_0 = U0*U0L_0 - 0.5d0*dexp(-lamda0*U0*U0)/dsqrt(pi*lamda0)
      U2L_0 = U0*U1L_0 + 0.5d0*U0L_0/lamda0
      U3L_0 = U0*U2L_0 + 1.d0*U1L_0/lamda0
      U4L_0 = U0*U3L_0 + 1.5d0*U2L_0/lamda0
      U5L_0 = U0*U4L_0 + 2.d0*U3L_0/lamda0
      U6L_0 = U0*U5L_0 + 2.5d0*U4L_0/lamda0

      gamma1_mul(1) = rho0
      gamma1_mul(2) = rho0*U1_0
      gamma1_mul(3) = rho0*V1_0
      gamma1_mul(4) = 0.5d0*rho0*(U2_0+V2_0+xi2)

      gamma2_mul(1) = a1barL*U1G_0+a2barL*U2G_0+a3barL*U1G_0*V1_0
     +               +a4barL*0.5d0*(U3G_0+U1G_0*V2_0+U1G_0*xi2)
     +               +a1barR*U1L_0+a2barR*U2L_0+a3barR*U1L_0*V1_0
     +               +a4barR*0.5d0*(U3L_0+U1L_0*V2_0+U1L_0*xi2)

      gamma2_mul(2) = a1barL*U2G_0+a2barL*U3G_0+a3barL*U2G_0*V1_0
     +               +a4barL*0.5d0*(U4G_0+U2G_0*V2_0+U2G_0*xi2)
     +               +a1barR*U2L_0+a2barR*U3L_0+a3barR*U2L_0*V1_0
     +               +a4barR*0.5d0*(U4L_0+U2L_0*V2_0+U2L_0*xi2) 

      gamma2_mul(3) = a1barL*U1G_0*V1_0+a2barL*U2G_0*V1_0
     +               +a3barL*U1G_0*V2_0
     +             +a4barL*0.5d0*(U3G_0*V1_0+U1G_0*V3_0+U1G_0*xi2*V1_0)
     +               +a1barR*U1L_0*V1_0+a2barR*U2L_0*V1_0
     +               +a3barR*U1L_0*V2_0
     +             +a4barR*0.5d0*(U3L_0*V1_0+U1L_0*V3_0+U1L_0*xi2*V1_0)

      gamma2_mul(4) = a1barL*0.5d0*(U3G_0+U1G_0*V2_0+U1G_0*xi2)
     +               +a2barL*0.5d0*(U4G_0+U2G_0*V2_0+U2G_0*xi2)
     +             +a3barL*0.5d0*(U3G_0*V1_0+U1G_0*V3_0+U1G_0*V1_0*xi2)
     +               +a4barL*0.25d0*(U5G_0+U1G_0*V4_0+U1G_0*xi4
     +             +2.d0*U3G_0*V2_0+2.d0*U3G_0*xi2+2.d0*U1G_0*V2_0*xi2)
     +               +a1barR*0.5d0*(U3L_0+U1L_0*V2_0+U1L_0*xi2)
     +               +a2barR*0.5d0*(U4L_0+U2L_0*V2_0+U2L_0*xi2)
     +             +a3barR*0.5d0*(U3L_0*V1_0+U1L_0*V3_0+U1L_0*V1_0*xi2)
     +               +a4barR*0.25d0*(U5L_0+U1L_0*V4_0+U1L_0*xi4
     +             +2.d0*U3L_0*V2_0+2.d0*U3L_0*xi2+2.d0*U1L_0*V2_0*xi2)


      gamma3_mul(1) = rhoL*U0GL+rhoR*U0LR
      gamma3_mul(2) = rhoL*U1GL+rhoR*U1LR
      gamma3_mul(3) = rhoL*U0GL*V1L+rhoR*U0LR*V1R
      gamma3_mul(4) = 0.5d0*(rhoL*(U2GL+U0GL*V2L+U0GL*xi2L)
     +                    +rhoR*(U2LR+U0LR*V2R+U0LR*xi2R)) 
   
      gamma4_mul(1) = a1L*U1GL+a2L*U2GL+a3L*U1GL*V1L
     +               +a4L*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +               +a1R*U1LR+a2R*U2LR+a3R*U1LR*V1R
     +               +a4R*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R)

      gamma4_mul(2) = a1L*U2GL+a2L*U3GL+a3L*U2GL*V1L
     +               +a4L*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +               +a1R*U2LR+a2R*U3LR+a3R*U2LR*V1R
     +               +a4R*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R)

      gamma4_mul(3) = a1L*U1GL*V1L+a2L*U2GL*V1L+a3L*U1GL*V2L
     +               +a4L*0.5d0*(U3GL*V1L+U1GL*V3L+U1GL*V1L*xi2L)
     +               +a1R*U1LR*V1R+a2R*U2LR*V1R+a3R*U1LR*V2R
     +               +a4R*0.5d0*(U3LR*V1R+U1LR*V3R+U1LR*V1R*xi2R)

      gamma4_mul(4) = a1L*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +               +a2L*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +               +a3L*0.5d0*(U3GL*V1L+U1GL*V3L+U1GL*V1L*xi2L)
     +               +a4L*0.25d0*(U5GL+U1GL*V4L+U1GL*xi4L
     +           +2.d0*U3GL*V2L+2.d0*U3GL*xi2L+2.d0*U1GL*V2L*xi2L)
     +               +a1R*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R)
     +               +a2R*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R)
     +               +a3R*0.5d0*(U3LR*V1R+U1LR*V3R+U1LR*V1R*xi2R)
     +               +a4R*0.25d0*(U5LR+U1LR*V4R+U1LR*xi4R
     +           +2.d0*U3LR*V2R+2.d0*U3LR*xi2R+2.d0*U1LR*V2R*xi2R)

      gamma5_mul(1) = gamma4_mul(1)
     +               +AL1*U0GL+AL2*U1GL+AL3*U0GL*V1L
     +               +AL4*0.5d0*(U2GL+U0GL*V2L+U0GL*xi2L) 
     +               +AR1*U0LR+AR2*U1LR+AR3*U0LR*V1R
     +               +AR4*0.5d0*(U2LR+U0LR*V2R+U0LR*xi2R)
      gamma5_mul(2) = gamma4_mul(2)
     +               +AL1*U1GL+AL2*U2GL+AL3*U1GL*V1L
     +               +AL4*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +               +AR1*U1LR+AR2*U2LR+AR3*U1LR*V1R
     +               +AR4*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R)

      gamma5_mul(3) = gamma4_mul(3)
     +               +AL1*U0GL*V1L+AL2*U1GL*V1L+AL3*U1GL*V1L
     +               +AL4*0.5d0*(U2GL*V1L+U0GL*V3L+U0GL*V1L*xi2L)
     +               +AR1*U0LR*V1R+AR2*U1LR*V1R+AR3*U1LR*V1R
     +               +AR4*0.5d0*(U2LR*V1R+U0LR*V3R+U0LR*V1R*xi2R)

      gamma5_mul(4) = gamma4_mul(4)
     +               +AL1*0.5d0*(U2GL+U0GL*V2L+U0GL*xi2L)
     +               +AL2*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +               +AL3*0.5d0*(U2GL*V1L+U0GL*V3L+U0GL*V1L*xi2L)
     +               +AL4*0.25d0*(U4GL+U0GL*V4L+U0GL*xi4L
     +               +2.d0*U2GL*V2L+2.d0*U2GL*xi2L+2.d0*U0GL*V2L*xi2L)
     +               +AR1*0.5d0*(U2LR+U0LR*V2R+U0LR*xi2R)
     +               +AR2*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R)
     +               +AR3*0.5d0*(U2LR*V1R+U0LR*V3R+U0LR*V1R*xi2R)
     +               +AR4*0.25d0*(U4LR+U0LR*V4R+U0LR*xi4R
     +             +2.d0*U2LR*V2R+2.d0*U2LR*xi2R+2.d0*U0LR*V2R*xi2R) 

      Ab1 = gamma1*gamma1_mul(1)+gamma2*gamma2_mul(1)
     +     +gamma3*gamma3_mul(1)+gamma4*gamma4_mul(1)
     +     +gamma5*gamma5_mul(1)

      Ab2 = gamma1*gamma1_mul(2)+gamma2*gamma2_mul(2)
     +     +gamma3*gamma3_mul(2)+gamma4*gamma4_mul(2)
     +     +gamma5*gamma5_mul(2)

      Ab3 = gamma1*gamma1_mul(3)+gamma2*gamma2_mul(3)
     +     +gamma3*gamma3_mul(3)+gamma4*gamma4_mul(3)
     +     +gamma5*gamma5_mul(3)

      Ab4 = gamma1*gamma1_mul(4)+gamma2*gamma2_mul(4)
     +     +gamma3*gamma3_mul(4)+gamma4*gamma4_mul(4)
     +     +gamma5*gamma5_mul(4)

      call invM(U0,V0,lamda0,cK,Ab1,Ab2,Ab3,Ab4,A1,A2,A3,A4)
ccccccc   Flux, F(1), F(2), F(3), F(4) cccccc 

      beta4 = tau*eta1
      beta5 = -tau*dt*eta+tau*beta4
      beta1 = dt-beta4
      beta2 = -tau*beta1+beta5
      beta3 = 0.5d0*dt*dt-dt*tau+tau*beta4

      beta1_mul(1) = rho0*U1_0      
      beta1_mul(2) = rho0*U2_0
      beta1_mul(3) = rho0*U1_0*V1_0
      beta1_mul(4) = 0.5d0*rho0*(U3_0+U1_0*V2_0+U1_0*xi2)

      beta2_mul(1) = gamma2_mul(2)
      beta2_mul(2) = a1barL*U3G_0+a2barL*U4G_0+a3barL*U3G_0*V1_0
     +              +a4barL*0.5d0*(U5G_0+U3G_0*V2_0+U3G_0*xi2)
     +              +a1barR*U3L_0+a2barR*U4L_0+a3barR*U3L_0*V1_0
     +              +a4barR*0.5d0*(U5L_0+U3L_0*V2_0+U3L_0*xi2)
      beta2_mul(3) = a1barL*U2G_0*V1_0+a2barL*U3G_0*V1_0
     +              +a3barL*U2G_0*V2_0
     +              +a4barL*0.5d0*(U4G_0*V1_0+U2G_0*V3_0+U2G_0*V1_0*xi2)
     +              +a1barR*U2L_0*V1_0+a2barR*U3L_0*V1_0
     +              +a3barR*U2L_0*V2_0
     +              +a4barR*0.5d0*(U4L_0*V1_0+U2L_0*V3_0+U2L_0*V1_0*xi2)

      beta2_mul(4) = a1barL*0.5d0*(U4G_0+U2G_0*V2_0+U2G_0*xi2)
     +              +a2barL*0.5d0*(U5G_0+U3G_0*V2_0+U3G_0*xi2)
     +              +a3barL*0.5d0*(U4G_0*V1_0+U2G_0*V3_0+U2G_0*V1_0*xi2)
     +              +a4barL*0.25d0*(U6G_0+U2G_0*V4_0+U2G_0*xi4+
     +              +2.d0*U4G_0*V2_0+2.d0*U4G_0*xi2+2.d0*U2G_0*V2_0*xi2)
     +              +a1barR*0.5d0*(U4L_0+U2L_0*V2_0+U2L_0*xi2)
     +              +a2barR*0.5d0*(U5L_0+U3L_0*V2_0+U3L_0*xi2)
     +              +a3barR*0.5d0*(U4L_0*V1_0+U2L_0*V3_0+U2L_0*V1_0*xi2)
     +              +a4barR*0.25d0*(U6L_0+U2L_0*V4_0+U2L_0*xi4+
     +             +2.d0*U4L_0*V2_0+2.d0*U4L_0*xi2+2.d0*U2L_0*V2_0*xi2) 

      beta3_mul(1) = A1*U1_0+A2*U2_0+A3*U1_0*V1_0
     +              +A4*0.5d0*(U3_0+U1_0*V2_0+U1_0*xi2)
      beta3_mul(2) = A1*U2_0+A2*U3_0+A3*U2_0*V1_0
     +              +A4*0.5d0*(U4_0+U2_0*V2_0+U2_0*xi2)
      beta3_mul(3) = A1*U1_0*V1_0+A2*U2_0*V1_0+A3*U1_0*V2_0
     +              +A4*0.5d0*(U3_0*V1_0+U1_0*V3_0+U1_0*V1_0*xi2)
      beta3_mul(4) = A1*0.5d0*(U3_0+U1_0*V2_0+U1_0*xi2)
     +              +A2*0.5d0*(U4_0+U2_0*V2_0+U2_0*xi2)
     +              +A3*0.5d0*(U3_0*V1_0+U1_0*V3_0+U1_0*V1_0*xi2)
     +              +A4*0.25d0*(U5_0+U1_0*V4_0+U1_0*xi4
     +              +2.d0*U3_0*V2_0+2.d0*U3_0*xi2+2.d0*U1_0*V2_0*xi2)


      beta4_mul(1) = -tau*(AL1*U1GL+AL2*U2GL+AL3*U1GL*V1L
     +                    +AL4*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +                    +AR1*U1LR+AR2*U2LR+AR3*U1LR*V1R
     +                    +AR4*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R))

      beta4_mul(2) = -tau*(AL1*U2GL+AL2*U3GL+AL3*U2GL*V1L
     +                    +AL4*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +                    +AR1*U2LR+AR2*U3LR+AR3*U2LR*V1R
     +                    +AR4*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R))

      beta4_mul(3) = -tau*(AL1*U1GL*V1L+AL2*U2GL*V1L+AL3*U1GL*V2L
     +                    +AL4*0.5d0*(U3GL*V1L+U1GL*V3L+U1GL*V1L*xi2L)
     +                    +AR1*U1LR*V1R+AR2*U2LR*V1R+AR3*U1LR*V2R
     +                    +AR4*0.5d0*(U3LR*V1R+U1LR*V3R+U1LR*V1R*xi2R))

      beta4_mul(4) = -tau*(AL1*0.5d0*(U3GL+U1GL*V2L+U1GL*xi2L)
     +                    +AL2*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +                    +AL3*0.5d0*(U3GL*V1L+U1GL*V3L+U1GL*V1L*xi2L)
     +                    +AL4*0.25d0*(U5GL+U1GL*V4L+U1GL*xi4L
     +                 +2.d0*U3GL*V2L+2.d0*U3GL*xi2L+2.d0*U1GL*V2L*xi2L)
     +                    +AR1*0.5d0*(U3LR+U1LR*V2R+U1LR*xi2R)
     +                    +AR2*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R)
     +                    +AR3*0.5d0*(U3LR*V1R+U1LR*V3R+U1LR*V1R*xi2R)
     +                    +AR4*0.25d0*(U5LR+U1LR*V4R+U1LR*xi4R
     +                +2.d0*U3LR*V2R+2.d0*U3LR*xi2R+2.d0*U1LR*V2R*xi2R))

      temp1_mul(1) = rhoL*U1GL+rhoR*U1LR
      temp1_mul(2) = rhoL*U2GL+rhoR*U2LR
      temp1_mul(3) = rhoL*U1GL*V1L+rhoR*U1LR*V1R
      temp1_mul(4) = 0.5d0*(rhoL*(U3GL+U1GL*V2L+U1GL*xi2L)
     +                   +rhoR*(U3LR+U1LR*V2R+U1LR*xi2R))

      temp2_mul(1) = a1L*U2GL+a2L*U3GL+a3L*U2GL*V1L
     +              +a4L*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +              +a1R*U2LR+a2R*U3LR+a3R*U2LR*V1R
     +              +a4R*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R)

      temp2_mul(2) = a1L*U3GL+a2L*U4GL+a3L*U3GL*V1L
     +              +a4L*0.5d0*(U5GL+U3GL*V2L+U3GL*xi2L)
     +              +a1R*U3LR+a2R*U4LR+a3R*U3LR*V1R
     +              +a4R*0.5d0*(U5LR+U3LR*V2R+U3LR*xi2R)

      temp2_mul(3) = a1L*U2GL*V1L+a2L*U3GL*V1L
     +              +a3L*U2GL*V2L
     +              +a4L*0.5d0*(U4GL*V1L+U2GL*V3L+U2GL*V1L*xi2L)
     +              +a1R*U2LR*V1R+a2R*U3LR*V1R
     +              +a3R*U2LR*V2R
     +              +a4R*0.5d0*(U4LR*V1R+U2LR*V3R+U2LR*V1R*xi2R)

      temp2_mul(4) = a1L*0.5d0*(U4GL+U2GL*V2L+U2GL*xi2L)
     +              +a2L*0.5d0*(U5GL+U3GL*V2L+U3GL*xi2L)
     +              +a3L*0.5d0*(U4GL*V1L+U2GL*V3L+U2GL*V1L*xi2L)
     +              +a4L*0.25d0*(U6GL+U2GL*V4L+U2GL*xi4L+
     +              +2.d0*U4GL*V2L+2.d0*U4GL*xi2L+2.d0*U2GL*V2L*xi2L)
     +              +a1R*0.5d0*(U4LR+U2LR*V2R+U2LR*xi2R)
     +              +a2R*0.5d0*(U5LR+U3LR*V2R+U3LR*xi2R)
     +              +a3R*0.5d0*(U4LR*V1R+U2LR*V3R+U2LR*V1R*xi2R)
     +              +a4R*0.25d0*(U6LR+U2LR*V4R+U2LR*xi4R+
     +              +2.d0*U4LR*V2R+2.d0*U4LR*xi2R+2.d0*U2LR*V2R*xi2R)

      F(1) = beta1*beta1_mul(1)+beta2*beta2_mul(1)+beta3*beta3_mul(1)
     +      +beta4*beta4_mul(1)+beta4*temp1_mul(1)
     +      -(beta5+tau*beta4)*temp2_mul(1)

      F(2) = beta1*beta1_mul(2)+beta2*beta2_mul(2)+beta3*beta3_mul(2)
     +      +beta4*beta4_mul(2)+beta4*temp1_mul(2)
     +      -(beta5+tau*beta4)*temp2_mul(2)

      F(3) = beta1*beta1_mul(3)+beta2*beta2_mul(3)+beta3*beta3_mul(3)
     +      +beta4*beta4_mul(3)+beta4*temp1_mul(3)
     +      -(beta5+tau*beta4)*temp2_mul(3)

      F(4) = beta1*beta1_mul(4)+beta2*beta2_mul(4)+beta3*beta3_mul(4)
     +      +beta4*beta4_mul(4)+beta4*temp1_mul(4)
     +      -(beta5+tau*beta4)*temp2_mul(4)

      FF(1) = F(1)/dt
      FF(2) = F(2)/dt
      FF(3) = F(3)/dt
      FF(4) = 0.d0!F(4)/dt
      FF(5) = F(4)/dt
c      write(*,*) F(1), F(2), F(3), F(4)  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
      return
      end
c    "A Gas-Kinetic BGK Scheme for the Navier-Stokes Equations
c    and Its connection with Artificial Dissipation and Godunov
c    Method", Kun Xu, Journal of Computational Physics, No 171, 
c    289-335 (2001), Appendix B

c    Hsiang-Hsu Wang Jan. 2007

      subroutine invM(U,V,lamda, K, b1, b2, b3, b4, a1, a2, a3, a4)
      real*8 U,V,lamda,K,b1,b2,b3,b4,a1,a2,a3,a4,R2,R3,R4

      R4 = 2.d0*b4-(U*U + V*V + 0.5d0*(K+2.d0)/lamda)*b1
      R3 = b3 - V*b1
      R2 = b2 - U*b1
      a4 = 4.d0*lamda*lamda/(K+2.d0)*(R4-2.d0*U*R2-2.d0*V*R3)
      a3 = 2.d0*lamda*R3-V*a4
      a2 = 2.d0*lamda*R2-U*a4
      a1 = b1-U*a2-V*a3-0.5d0*a4*(U*U+V*V+0.5d0*(K+2.d0)/lamda)


c      write(*,*) U, V, lamda0, K
      return
      end
