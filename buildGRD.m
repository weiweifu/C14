function [GRD] = buildGRD(F,varargin)
% function [GRD] = buildGRD(F,varargin)
% returns a structure, MET, containing geometry for the configuration 
% Inputs:  ENCODED HEREIN:
%  resolution is used to construct filename of the netcdf output file 
%      resolution can be '240km', '120km', '60km', '30km', and '15km' 
% 
%  configuration dimensions: ny,nz
%
%  marginal_mask is an (additional) mask for the marginal seas to
%      be excluded from the ocean configuration for the operators to be 
%      used.
%  Notes: 1) Dimensions output are in meters.
%         2) Variables are generally output in 3_D, even if this means they are
%            replicated. Exceptions are 
%              KMT, which is 2-D, and
%              iocn, which is a vector of 1-D indices into the ocean cells
%              of the 3-D grid geometry.
%         3) Refer to the POP Reference Manual, Chapter 3, for figures
%            and further text depicting the geometry definitions.
%            In particular, see figures 3.1, 3.2, and 3.3. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%  CONFIGURATION DEFINITIONS - to be modified for a different configuration:
  test_dir     = F.IRFdir;
  % test_case    = ['worldOcean_QU_' resolution '/restarts/'];
  % ref_dir      = [test_dir test_case ];
  % D            = dir([ref_dir 'restart*nc']);
  % ref_file     = D.name;
  % The netcdf file from which to get the geometry
  ref_dir    = test_dir;
  ref_file   = F.IRFfile;

  % The directory to put the GRD struct
  mask_dir = F.irf_mask_dir;

  %
  % read in the metric coefficients frpm the netcdf file
  ncid = netcdf.open([ref_dir,ref_file],'NOWRITE');

  i1 = netcdf.inqDimID(ncid,'nEdges');       % nEdges
  i2 = netcdf.inqDimID(ncid,'maxEdges');     % maxEdges
  i3 = netcdf.inqDimID(ncid,'nCells');       % nCells
  i4 = netcdf.inqDimID(ncid,'nVertLevels');  % nVertLevels

  [dimname, nuy] = netcdf.inqDim(ncid,i1);  % nEdges
  [dimname, muy] = netcdf.inqDim(ncid,i2);  % maxEdges
  [dimname, ny]  = netcdf.inqDim(ncid,i3);  % nCells
  [dimname, nz]  = netcdf.inqDim(ncid,i4);  % nVertLevels
  nzz            = nz + 1;

  MASK_MARGINAL = ones(ny,nzz);       % include all except marginal seas

  
  % use the KMT (index of the bottom wet cell (land at the surface = 0)
  % to create the land sea mask  
  varid = netcdf.inqVarID(ncid,'maxLevelCell'); % number of wet cells above seafloor
  KMT   = netcdf.getVar(ncid,varid,'double');

  if length(KMT) ~= ny || max(KMT) ~= nz; 
    error('Please check the dimensions!')
  end
%  nz = max(KMT(:)); nzz = nz+1;        % Offline variables have a seafloor;
                                        % k index increased by 1;
  fprintf('Loading the landmask and metric coefficients..\n');
  MASK = zeros(ny,nzz);
  for k = 1:nzz
    MASK(:,k) = (k <= KMT(:));
  end % all k's 
  
  varid = netcdf.inqVarID(ncid,'temperature'); % number of wet cells above seafloor
  TEMP  = netcdf.getVar(ncid,varid,'double');
  % check expected land-sea masks 
  [nyck,nzzck] = size(MASK);
  if (nyck ~=ny) || (nzzck ~=nzz)
    display('dimensions not as expected in buildMET');
    keyboard;
  elseif ~isempty(find(TEMP'.*MASK(:,1:nz)<=-1e34))
    display('A dry point is set to wet!!');
  elseif numel(find(MASK(:,1:nz)==0)) ~= numel(find(TEMP==-1e34))
    display('dry points in MASK is not equal to those in temperature!!')
  else
    display('MASK is correct...');
  end
  
  %  mask out the small and irrelevant closed basins
  % MASK = MASK .* MASK_MARGINAL;
  % get reference index for Cell, Edge and Vertex, which is necessary to reference
  % neighbours in the Voronoi regions and Delaunay triangles

  % cellsOnCell List of cells that neighbor each cell.
  varid = netcdf.inqVarID(ncid,'cellsOnCell');
  cOnC  = netcdf.getVar(ncid,varid,'double' );
  % edgesOnCell List of edges that border each cell.
  varid = netcdf.inqVarID(ncid,'edgesOnCell');
  eOnC  = netcdf.getVar(ncid,varid,'double' );
  % nEdgesOnCell Number of edges that border each cell. 
  varid = netcdf.inqVarID(ncid,'nEdgesOnCell');
  neOnC  = netcdf.getVar(ncid,varid,'double' );
  % verticesOnCell List of vertices that border each cell.
  varid = netcdf.inqVarID(ncid,'verticesOnCell');
  vOnC  = netcdf.getVar(ncid,varid,'double' );

  % nEdgesOnEdge Number of edges that surround each of the cells that straddle each edge.
  % These edges are used to reconstruct the tangential velocities.
  % cellsOnEdge List of cells that straddle each edge.
  % verticesOnEdge List of vertices that straddle each edge.
  % edgesOnEdge List of edges that border each of the cells that straddle each edge.
  % weightsOnEdge Reconstruction weights associated with each of the edgesOnEdge.

  % edgesOnVertex List of edges that share a vertex as an endpoint.
  % cellsOnVertex List of cells that share a vertex.
  % 
  % get horizontal length of Edges measuring the distance between nCells   
  % dcEdge Length of each edge, computed as the distance between cellsOnEdge.
  varid = netcdf.inqVarID(ncid,'dcEdge');
  dcEdge   = netcdf.getVar(ncid,varid,'double' );
  dcEdge   = dcEdge(:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %
  % get horizontal length of Edges measuring the distance between nVertices   
  % dvEdge Length of each edge, computed as the distance between verticesOnEdge. 
  varid = netcdf.inqVarID(ncid,'dvEdge');
  dvEdge   = netcdf.getVar(ncid,varid,'double' );
  dvEdge   = dvEdge(:,ones(nzz,1));  % extend to 3_D, with added layer at the bottom
  %  
  % get layerThickness for resting state, vertical thickness of the nCells (same in all horizontal
  % directions)
  varid = netcdf.inqVarID(ncid,'restingThickness');
  DZT(:,:) = netcdf.getVar(ncid,varid,'double' ); % unit [m]
  DZT(nzz,:) = DZT(nz,:);            % small non-zero depth for added layer
  DZT = permute(DZT,[2 1]);
  %  
  % get refBottomDepth: Reference depth of ocean for each vertical level. 
  % Used in ’z-level’ type runs.
  varid = netcdf.inqVarID(ncid,'refBottomDepth');
  refBottomDepth  = netcdf.getVar(ncid,varid,'double' ); % unit [m]
  ZT    = netcdf.getVar(ncid,varid,'double' ); % unit [m]
  ZT(nzz) = 5600;
  ZT = repmat(ZT',[ny 1]);    % extend to 3_D
  %
  % lat and long, at the middle of the nCells
  varid   = netcdf.inqVarID(ncid,'lonCell');
  lonCell = netcdf.getVar(ncid,varid,'double' );
  varid   = netcdf.inqVarID(ncid,'latCell');
  latCell = netcdf.getVar(ncid,varid,'double' );
  lonCell3D = repmat(lonCell,[1 nzz]);   % extend to 3_D
  latCell3D = repmat(latCell,[1 nzz]);   % extend to 3_D
  %  
  % AREA of Cells, area of each Cell
  % areaCell Area of each cell in the primary grid.
  varid    = netcdf.inqVarID(ncid,'areaCell');
  areaCell = netcdf.getVar(ncid,varid,'double' );
  areaCell = repmat(areaCell,[1 nzz]);   % extend to 3_D
  %  
  % volume of the nCell (m^3)
  volCell = areaCell.*DZT;
  %
  % lat and long, for Edges
  varid   = netcdf.inqVarID(ncid,'lonEdge');
  lonEdge = netcdf.getVar(ncid,varid,'double' );
  varid   = netcdf.inqVarID(ncid,'latEdge');
  latEdge = netcdf.getVar(ncid,varid,'double' );
  % lat and long, for Vertices
  varid     = netcdf.inqVarID(ncid,'lonVertex');
  lonVertex = netcdf.getVar(ncid,varid,'double' );
  varid     = netcdf.inqVarID(ncid,'latVertex');
  latVertex = netcdf.getVar(ncid,varid,'double' );
  %  
  % AREA of Triangles, area of each Triangle
  % areaTriangle Area of each cell (triangle) in the dual grid.
  varid        = netcdf.inqVarID(ncid,'areaTriangle');
  areaTriangle = netcdf.getVar(ncid,varid,'double' );
   
  %
  netcdf.close(ncid);
  
  % create the structure 
  GRD.cOnC         = cOnC;
  GRD.neOnC        = neOnC;
  GRD.lonCell      = lonCell;
  GRD.latCell      = latCell;
  GRD.areaCell     = areaCell;
  GRD.lonCell3D    = lonCell3D;
  GRD.latCell3D    = latCell3D;
  % GRD.areaTriangle = areaTriangle;

  GRD.MASK = MASK;
  GRD.iocn = find(GRD.MASK(:)~=0); % 1-D indices into the 3-D grid structure
                                   % for ocean T_cells.
  GRD.KMT  = KMT;
  GRD.TEMP = TEMP;
  GRD.ZT   = ZT;
  GRD.DZT  = DZT;
  GRD.ny   = ny;
  GRD.nz   = nz;
  GRD.nzz  = nz+1;
  GRD.volCell = volCell;
  GRD.refBottomDepth = refBottomDepth;

  GRD.grain = 1;                   % full sized MET
  if nargin >1
    disp(['saving GRD to ' mask_dir]);
    save(sprintf('%sMSK_%s.mat',mask_dir,F.resolution),'GRD','-v7.3');
  end
  return 

end % function buildGRD
