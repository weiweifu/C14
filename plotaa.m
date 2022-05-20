function plotaa(GRD,dd2,nlev,method)
% function plotaa(GRD,dd2,nlev,method)
% method can be one of the following:
%  'nearest'   - Nearest neighbor interpolation
%  'linear'    - Linear interpolation (default)
%  'natural'   - Natural neighbor interpolation
%  'cubic'     - Cubic interpolation (2D only)
%  'v4'        - MATLAB 4 griddata method (2D only)
 
addpath ~/MatUtils/misc/; 
% [vq,xq,yq] = interp2xy(GRD.lonCell*180/pi,GRD.latCell*180/pi,GRD.TEMP(10,:));
if nargin > 3
  [vq,xq,yq] = interp2xy(GRD.lonCell*180/pi,GRD.latCell*180/pi,real(dd2(:,nlev)),method);
else
  [vq,xq,yq] = interp2xy(GRD.lonCell*180/pi,GRD.latCell*180/pi,real(dd2(:,nlev)));
end
plotxy('proj','miller','dlon',xq,'dlat',yq,'data',vq,'plt','pcolor');
