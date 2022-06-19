
clear all;
ibeg = 1;
iend = 100;
jbeg = 1;
jend = 100;
ibuf = 2;
jbuf = 2;
Nx   = iend+ibuf*2.0;
Ny   = jend+jbuf*2.0;

rho(Nx,Ny)=0.0;
fid = fopen('rho0078.dat','r');
rhotemp = fscanf(fid,'%g',Nx*Ny);
k=1;
for i=1:1:Nx
	for j=1:1:Ny
		rho(i,j)=(rhotemp(k));
		k=k+1;
	end
end
rho=transpose(rho);
mesh(rho,'FaceColor','interp')
view(2)
