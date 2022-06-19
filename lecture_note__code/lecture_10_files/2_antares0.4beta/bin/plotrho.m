%function showrho(den)
% Load Data
dataload
%
   den =-3
   figure(1)
   for i=1:1:Nx
   for j=1:1:Ny
   for k=1:1:Nz
   rho2(j,i,k)=log10(rho(i,j,k));
   end
   end
   end
%   slice(rho2,0,(Nx+2),3)
   slice(rho2,0,(Nx+2),34)
   shading flat
   hold on
   isosurface(rho2,den)
   isocaps(rho2,den)
   camlight
   camlight right
   lighting phong
   daspect([1,1,1])
   set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
   set(gca,'Color','black','XColor','white','YColor','white','ZColor','white')
   box on
   view(50, 30)
%   caxis([-3.5 -0.5]);
%   mesh(xx,yy,rho,'FaceColor','interp');
   colorbar;
%   axis(vxs);
   axis([0 Nx 0 Ny 0 Nz]);
   axis equal;
   title('Density');
   hold off
