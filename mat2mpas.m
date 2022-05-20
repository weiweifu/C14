function [out] = mat2mpas(in,varargin)
% function [out] = mat2mpas(in,varargin)
% convert 1-D vector to mpas tracer matrix (nVertLevels, nCells)
% to check variables. Note the last layer (nzz) is all 0's for tracers. 
% if GRD is an argument, convert in(iocn) to (nzz,nCells)
 
 if nargin > 1
   if numel(in)~= numel(varargin{1}.iocn)
     error('check input matrix');
   end
   [ny,nzz] = size(varargin{1}.MASK);
   out = repmat(nan,[ny,nzz-1]);
   out(varargin{1}.iocn) = in;
 else
   ik = size(in,1);
   jk = ik/kN;
 
   out = reshape(in,[jk kN]);
   out(:,kN) = [];
 
   out = permute(out,[2 1]);
 end
 end
 
