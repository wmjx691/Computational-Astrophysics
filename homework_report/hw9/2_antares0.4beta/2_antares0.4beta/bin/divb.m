
% Load Data
dataload
%

 dbxdx = ([bx(:,3) bx(:,3:Nx) bx(:,Nx)]-[bx(:,1) bx(:,1:Nx-2) bx(:,Nx-2)])/(2.0*dx);
 dbydy = ([by(3,:);by(3:Ny,:);by(Ny,:)]-[by(1,:);by(1:Ny-2,:);by(Ny-2,:)])/(2.0*dy);
 dbydy = ([by(:,3) by(:,3:Ny) by(:,Ny)]-[by(:,1) by(:,1:Ny-2) by(:,Ny-2)])/(2.0*dy);
 dbxdx = ([bx(3,:);bx(3:Nx,:);bx(Nx,:)]-[bx(1,:);bx(1:Nx-2,:);bx(Nx-2,:)])/(2.0*dx);

 divB = dbxdx + dbydy;

   mesh(xx,yy,divB,'FaceColor','interp');
   colorbar;
   axis(vxs);
