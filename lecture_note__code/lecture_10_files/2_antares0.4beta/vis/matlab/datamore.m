%================================================================
%
% Evaluate more data information
%
% Kuo-Chuan Pan Feb 29,2008
%---------------------------------------------------------------
% load data
  dataload;
% Velocity
  vx=px./rho;
  vy=py./rho;
% calculate r and r^2  
  rsq = x.^2 +y.^2;
  r   = sqrt(rsq);
% calculate velocity in Polar coordinate 
  vr  =( vx.*x + vy.*y)./r; 
  vs  =(-vx.*y + vy.*x)./r; 


