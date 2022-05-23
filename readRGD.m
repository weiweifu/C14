% the file containing the output of the regridding calculation
fname = 'W30to60_240km.nc';
fname = 'W30to60_240km_181106.nc';
fname = 'W240km_30to60_181106.nc';
% src_grid_rank,dst_grid_rank
% The number of dimensions of the source grid. Currently this can only be 1 or 2.
% Where 1 indicates an unstructured grid and 2 indicates a 2D logically rectangular grid.

% source grid information
% n_a: The number of source cells.
% nv_a: The maximum number of corners (i.e. vertices) around a source cell. 
% If a cell has less than the maximum number of corners, then the remaining 
% corner coordinates are repeats of the last valid corner's coordinates.

% The longitude coordinates of the centers of each source cell.
rgd.xc_a = ncread(fname,'xc_a',[1],[inf]);
% The latitude coordinates of the centers of each source cell.
rgd.yc_a = ncread(fname,'yc_a',[1],[inf]);
% The longitude coordinates of the corners of each source cell.
rgd.xv_a = ncread(fname,'xv_a',[1 1],[inf inf]);
% The latitude coordinates of the corners of each source cell.
rgd.yv_a = ncread(fname,'yv_a',[1 1],[inf inf]);
% The mask for each source cell. A value of 0, indicates that the cell is masked.
rgd.mask_a = ncread(fname,'mask_a',[1],[inf]);
% The area of each source cell. This quantity is either from the source 
% grid file or calculated by ESMF_RegridWeightGen. When a 
% non-conservative regridding method (e.g. bilinear) is used, the area is set to 0.0.
rgd.area_a = ncread(fname,'area_a',[1],[inf]);
% The number of cells along each dimension of the source grid. 
% For unstructured grids this is equal to the number of cells in the grid.
rgd.src_grid_dims = ncread(fname,'src_grid_dims',[1],[inf]);
rgd.n_a = size(rgd.xv_a,2);
rgd.nv_a = size(rgd.xv_a,1);

% The longitude coordinates of the centers of each destination cell.
rgd.xc_b = ncread(fname,'xc_b',[1],[inf]);
% The latitude coordinates of the centers of each destination cell.
rgd.yc_b = ncread(fname,'yc_b',[1],[inf]);
% The longitude coordinates of the corners of each destination cell.
rgd.xv_b = ncread(fname,'xv_b',[1 1],[inf inf]);
% The latitude coordinates of the corners of each destination cell.
rgd.yv_b = ncread(fname,'yv_b',[1 1],[inf inf]);
% The mask for each destination cell. A value of 0, indicates that the cell is masked.
rgd.mask_b = ncread(fname,'mask_b',[1],[inf]);
% The area of each destination cell. This quantity is either from the destination 
% grid file or calculated by ESMF_RegridWeightGen. When a 
% non-conservative regridding method (e.g. bilinear) is used, the area is set to 0.0.
rgd.area_b = ncread(fname,'area_b',[1],[inf]);
% The number of cells along each dimension of the destination grid. 
% For unstructured grids this is equal to the number of cells in the grid.
rgd.dst_grid_dims = ncread(fname,'dst_grid_dims',[1],[inf]);
rgd.n_b = size(rgd.xv_b,2);
rgd.nv_b = size(rgd.xv_b,1);

% The position in the source grid for each entry in the regridding matrix.
rgd.col = ncread(fname,'col',[1],[inf]);
% The position in the destination grid for each entry in the weight matrix.
rgd.row = ncread(fname,'row',[1],[inf]);
% The weight for each entry in the regridding matrix.
rgd.S = ncread(fname,'S',[1],[inf]);
%  the fraction of each source cell that participated in the regridding
rgd.frac_a = ncread(fname,'frac_a',[1],[inf]);
% the fraction of each destination cell that participated in the regridding.
rgd.frac_b = ncread(fname,'frac_b',[1],[inf]);
% The number of entries in the regridding matrix
rgd.n_s = size(rgd.S,1);
