 function [out] = mpas2mat(in,varargin)
 % function [out] = mpas2mat(in,varargin)
 % convert mpas tracer matrix (nVertLevels, nCells) to matlab
 % multiplying form by adding 1 layer of 0's to the bottom
 % if iocn is an argument, only wet points are considered.
 
   ik = size(in,1);
   jk = size(in,2);
   % disp(['Input shape: ' mat2str(size(in))]);

   temp = in;

   temp = permute(temp,[2 1]);
   % disp(['Shape after permute: ' num2str(size(temp))]);
 
   % temp(ik+1,:) = zeros(1,jk);
   temp(:,ik+1) = zeros(1,jk);
   % disp(['Input shape now:      ' num2str(size(temp))]);
 
   if nargin >1
     out = temp(varargin{1}.iocn);
   else
     out = temp(:);
   end
 end
