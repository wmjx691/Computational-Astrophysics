      integer ibeg,iend,ibuf
      integer jbeg,jend,jbuf
      integer kbeg,kend,kbuf
      integer i,j,k,nprnt
      integer istep,nstop
      integer nx,ny,nz,nd
      real*8 x,xl,y,yl,z,zl,r
      real*8 dx,dy,dz,dt,t,tf
      real*8 gam,dK,CV,rho,vx,vy,vz,en,iso_snd
      real*8 xmin,xmax,ymin,ymax,zmin,zmax,cfl
      real*8 fio_ppm,fio_bin_i,fio_bin_o
      real*8 pt
      real*8 smalld,smalle,smallt
      logical d1d,d2d,d3d,mhd,isothermal
      logical cartesian, spherical, cylindrical
      common /geom/ x,xl,y,yl,z,zl
  
