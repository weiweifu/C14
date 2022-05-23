function indx = fd(i,GRD)
% function indx = fd(i,GRD)
% find the neighbor of grid (i), a vector is returned
  indx = GRD.cOnC(1:GRD.neOnC(i),i);
  % indx = nonzeros(indx);
  indx = setdiff(indx,[0]);
  indx = [i;indx];
end
