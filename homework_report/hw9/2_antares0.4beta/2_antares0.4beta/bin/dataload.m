%===================================================
%
%  Load Data
%
%  Kuo-Chuan Pan Feb 29,2008
%
%--------------------------------------------------
clear all;

ibeg = 1;
iend = 60;
ibuf = 2;
jbeg = 1;
jend = 60;
jbuf = 2;
kbeg = 1;
kend = 60;
kbuf = 2;
Nx   = iend+ibuf*2.0;
Ny   = jend+jbuf*2.0;
Nz   = kend+kbuf*2.0;

xmin = -50;
xmax =  50;
ymin = -50;
ymax =  50;
zmin = -50;
zmax =  50;

dx  = (xmax-xmin)/(Nx-2*ibuf);
dy  = (ymax-ymin)/(Ny-2*jbuf);
dz  = (zmax-zmin)/(Nz-2*kbuf);

 xl  = xmin-dx:dx:xmax+2*dx;
 yl  = ymin-dy:dy:ymax+2*dy;
 zl  = zmin-dz:dz:zmax+2*dz;
 x  = xmin-1.5*dx:dx:xmax+1.5*dx;
 y  = ymin-1.5*dy:dy:ymax+1.5*dy;
 z  = zmin-1.5*dz:dz:zmax+1.5*dz;

[yy,xx]=meshgrid(y,x);
	 
vxs  = [xmin,xmax,ymin,ymax];


rho(Nx,Ny,Nz)=0.0;
 px(Nx,Ny,Nz)=0.0;
 py(Nx,Ny,Nz)=0.0;
 pz(Nx,Ny,Nz)=0.0;
ene(Nx,Ny,Nz)=0.0;
 bx(Nx,Ny,Nz)=0.0;
 by(Nx,Ny,Nz)=0.0;
 bz(Nx,Ny,Nz)=0.0;

fid = fopen('den0042.dat','r');
rhotemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('px0042.dat','r');
pxtemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('py0042.dat','r');
pytemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('pz0042.dat','r');
pytemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('ene0042.dat','r');
enetemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('bx0042.dat','r');
bxtemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('by0042.dat','r');
bytemp = fscanf(fid,'%g',Nx*Ny*Nz);

fid = fopen('bz0042.dat','r');
bztemp = fscanf(fid,'%g',Nx*Ny*Nz);

n=1;
for i=1:1:Nx
	for j=1:1:Ny
		for k=1:1:Nz
		rho(i,j,k)=(rhotemp(n));
		 px(i,j,k)= (pxtemp(n));
		 py(i,j,k)= (pytemp(n));
		 pz(i,j,k)= (pytemp(n));
		ene(i,j,k)=(enetemp(n));
		 bx(i,j,k)= (bxtemp(n));
		 by(i,j,k)= (bytemp(n));
		 bz(i,j,k)= (bztemp(n));
		n=n+1;
	        end
	end
end
clear rhotemp;
clear  pxtemp;
clear  pytemp;
clear  pztemp;
clear enetemp;
clear  bxtemp;
clear  bytemp;
clear  bztemp;

%rho=transpose(rho);
% px=transpose(px);
% py=transpose(py);
%ene=transpose(ene);
% bx=transpose(bx);
% by=transpose(by);
% bz=transpose(bz);

%figure(1)
%mesh(rho,'FaceColor','interp')
%view(2)

%figure(2)
%mesh(px,'FaceColor','interp')
%view(2)








