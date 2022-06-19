%function showrho(den)
% Load Data
dataload
%
   den =0.5
   figure(1)
%   rho2=log10(rho);
   rho2=rho(:,:,30);
%   slice(rho2,0,(Nx+2),3)
%   shading flat
%   hold on
%   isosurface(rho2,den)
%   isocaps(rho2,den)
%   camlight
%   camlight right
%   lighting phong
%   daspect([1,1,1])
%   set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
%   set(gca,'Color','black','XColor','white','YColor','white','ZColor','white')
%   box on
%   view(50, 30)
%   caxis([-3.5 -0.5]);
   mesh(xx,yy,rho2,'FaceColor','interp');
   colorbar;
%   axis(vxs);
   axis([0 Nx 0 Ny 0 Nz]);
   axis equal;
   title('Density');
   hold off
