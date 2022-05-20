function [vq,varargout] = interp2xy(x,y,v,method)
%  'nearest'   - Nearest neighbor interpolation
%  'linear'    - Linear interpolation (default)
%  'natural'   - Natural neighbor interpolation
%  'cubic'     - Cubic interpolation (2D only)
%  'v4'        - MATLAB 4 griddata method (2D only)
 nout = max(nargout,1)-1;
% function vq = interp2latlon(x,y,v,xq,yq)
x1 = 0.5:1:359.5;
y1 = -89.5:89.5;
[xq,yq] = meshgrid(x1,y1);
s = {xq,yq};
for i=1:nout
  varargout{i} = s{i}; 
end
% dt = delaunayTriangulation(X,Y);
% F = griddedInterpolant(squeeze(x),squeeze(y),squeeze(v));
if nargin > 3
  disp(['using method ' method]);
else
  method = 'natural';
end
vq = griddata(x,y,v,xq,yq,method);
% vq = F(xq,yq);
