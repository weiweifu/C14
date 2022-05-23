function indx = fdd(i,GRD)
% find the neighbor of grid (i), a vector is returned
% then find the neighbor of the neighbor of grid (i), a vector is returned
  idx0 = fd(i,GRD);
  indx = [i];
  % disp(mat2str(idx0));
  % indx is the neighbors of i (including i itself)
  for j=1:length(idx0)
    % loop over the neighbors of i and find their neighbors
    k = idx0(j);
    % disp(['now find neighbors for ' int2str(k)]);
    tmp = fd(k,GRD);
    % disp([num2str(k) ' : ' mat2str(tmp)]);
    indx = union(indx,tmp);
  end

% function indx = fd(i,GRD)
% % find the neighbor of grid (i), a vector is returned
%   indx = GRD.cOnC(1:GRD.neOnC(i),i);
%   indx = nonzeros(indx);
