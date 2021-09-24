function [LUMP,SPRAY,SL, MET_coarse] = lumpX(MET) 
%  create coarse graining volume-weighted operators
%
%  B = LUMP*A*SPRAY
%  use the unix pipeline  philosophy 
%  read from right to left: 
%  (1) spray the coarse grained field into a high res field
%  (2) apply the differential operator on the high res field
%  (3) lump the high res field back into the coarse grained field
%
%  Returns the LUMP and SPRAY operators, and a compatible version of MET. 
%  The LUMP and SPRAY operators returned are to be applied to iocn-sized
%  operators or vectors.
%  SL is a sum operator, not volume weighted, and is applied to 3_D 
%  full-dimensioned objects, like volume, where the sum, rather than the
%  average is wanted.
%
%  To further condense operators, go through this procedure multiple
%  times. 
%
% INPUTS:
%    MET                     % geometric data for the size of grid to be 
%                            % condensed (full size -> condensed-by-2 (X2);
%                            % X2 -> X4)
% OUTPUTS:
%    LUMP                    % LUMP operator, to be applied to iocn-sized
%                            % advection and diffusion operators, or to 
%                            % iocn-sized data vectors
%    SPRAY                   % SPRAY operator, to be appled in conjunction
%                            % with LUMP for iocn_sized advection and 
%                            % diffusion operators.
%    SL                      % summing operator, to be used on full
%                            % grid-sized data, where a sum rather than
%                            % an average is needed (as for volume).
%                            % SL is not volume weighted.
%    MET_coarse              % geometry for the coarse-grained resulting
%                            % operators. The fields present are:
%      iocn                  % vector of 1-D indices of ocean cells in MASK
%      MASK                  % ocean IFF any ocean cells in the lumping
%      VOL                   % sum of the ocn cells lumped together
%      DZT                   % same as in MET (all grid cells in one layer
%                            % have the same dimension
%      DZU                   % same as in MET
%      ZT                    % same as in MET
%      ZU                    % same as in MET
%      TAREA                 % sum of the ocn cells lumped together
%      TLAT                  % latitude of the bottom left grid-cell being
%                            % lumped
%      TLONG                 % longitude of the bottom left grid-cell being
%                            % lumped
%      KMT                   % re-made from the coarse-grained MASK
%      grain                 % 2*MET.grain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  dVt = MET.VOL;
  msk = MET.MASK;
  iocn = find(msk(:)==1);
  ilnd = find(msk(:)==0);
  
  % build an operator that sums 4x4 horizontal groups of grid boxes
  [ny,nx,nz] = size(dVt);
  ex = ones(nx,1); Ix = d0(ex);
  ey = ones(ny,1); Iy = d0(ey);
  ez = ones(nz,1); Iz = d0(ez);
  
  % sum every other grid box in the x direction
  ex2 = ones(nx/2,1); Ax = kron(d0(ex2),[1,1]);
  % sum every other grid box in the y direction
  ey2 = ones(ny/2,1); Ay = kron(d0(ey2),[1 1]);
  
  % Left multiply sums 4x4 group of grid boxes
  SL = kron(Iz,kron(Ax,Ay)); % transpose sprays lump into 4x4 group of grid boxes
  
  % compute the volumes of the coarse grained grid boxes
  dVt(ilnd) = 0; % zero out the volumes of the land points
  dV = d0(dVt(:)); 
  %
  % use the unix pipeline philosophy
  % read from right to left:
  % (1) multiply by the grid box volume
  % (2) sum 4x4 group
  % (3) divide by the total volume of the 4x4 group 
  LUMP_full = d0(1./(SL*dVt(:)))*SL*dV;          % volume-weighted averages
  %
  SPRAY_full = SL'; % transpose of the sum operator is the spray operator
  %  
  % create a coarse grained land-ocean mask
  mask_coarse = zeros(ny/2,nx/2,nz);
  mask_coarse(:) = (LUMP_full*(msk(:))~=0);     
  iocn_coarse = find(mask_coarse(:)==1);
  MET_coarse.iocn = iocn_coarse;
  
  % lumped volume matrix ( sum of the volumes of the fine-grained input )
  dVt_coarse = 0*mask_coarse;
  dVt_coarse(:) = SL*dVt(:); 
 
  % keep only the ocn points on both the fine and the coarse grid
  SPRAY = SPRAY_full(:,iocn_coarse);
  SPRAY = SPRAY(iocn,:);
  LUMP = LUMP_full(:,iocn);
  LUMP = LUMP(iocn_coarse,:);
  
  MET_coarse.MASK = mask_coarse;
  MET_coarse.VOL = dVt_coarse;
  
  [ny_c,nx_c,nz_c] = size(mask_coarse);
  
  % DZT and DZU are the same at all horizontal locations, as are
  % ZT, depth of T_points (positive, from the surface, mid-T_cell), and
  % ZU, depth of the top of T_cells (positive, from the surface)
  MET_coarse.DZT = MET.DZT(1:ny_c,1:nx_c,1:nz_c);
  MET_coarse.DZU = MET.DZU(1:ny_c,1:nx_c,1:nz_c);
  MET_coarse.ZT  = MET.ZT(1:ny_c,1:nx_c,1:nz_c);
  MET_coarse.ZU  = MET.ZU(1:ny_c,1:nx_c,1:nz_c);
  MET_coarse.DXT = MET.DXT(1:ny_c,1:nx_c,1:nz_c);
  MET_coarse.DYT = MET.DYT(1:ny_c,1:nx_c,1:nz_c);
  
  % area is the sum of the areas of the ocean cells
  AREA = MET.TAREA;
  AREA(ilnd) = 0;                      % set land cell area = 0 for summing
  AREA_coarse = 0*mask_coarse;
  AREA_coarse(:) = SL*AREA(:);
  MET_coarse.TAREA = AREA_coarse;
 
  % TLAT and TLONG are picked components of the fine-grained grid
  % (Hint: averaging does not give desireable results because of the
  % bends in longitude in the far North.)
  %
  % pick every other grid box in the x direction
  Axp = kron(d0(ex2),[1,0]);
  % pick every other grid box in the y direction
  Ayp = kron(d0(ey2),[1 0]);
  
  % Left multiply sums 4x4 group of grid boxes
  pickL = kron(Iz,kron(Axp,Ayp));      % pick the value in the first of 
                                       % the grid boxes  
  TLAT_coarse      = 0*mask_coarse;
  TLAT_coarse(:)   = (pickL*MET.TLAT(:));
  MET_coarse.TLAT  = TLAT_coarse;
 
  TLONG_coarse     = 0*mask_coarse;
  TLONG_coarse(:)  = (pickL*MET.TLONG(:));
  MET_coarse.TLONG = TLONG_coarse;
 
  % make KMT from the coarse mask
  KMT = zeros(ny_c,nx_c) + nz_c;
  for k = 0:(nz_c-1)
      land = find(MET_coarse.MASK(:,:,(k+1) ) == 0 );
      KMT(land) = min(KMT(land),k);    % set level if not set at highter level
  end % for each layer above the bott
  MET_coarse.KMT = KMT;
 
  MET_coarse.grain = MET.grain*2;
  
  return
end % function lumpX

