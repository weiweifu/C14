% set up transport and grid file
rdir = '/DFS-L/DATA/moore/weiweif/MPAS-BGC/';
%res = '240km';
res = '60km';
% mach = 'cori';
mach = 'gp';
model = 'iage';

TTM = {'A','H','D','T'};
fnp = '_iocn_mo_00_';

switch lower(mach)
  case 'cori'
    dir0 = '/global/cfs/cdirs/e3sm/weiweif/MPAS_IRF/Mask_TTM/';
    Tname = [dir0 'T_iocn_mo_00_' res '.mat'];
    GRDname = ['MSK_' res '.mat'];
  case 'gp'
    % dir0 = '/DFS-L/DATA/moore/weiweif/MPAS-BGC/Data/';
    dir0 = [rdir 'Mask_TTM/'];
    % Tname = ['T_iocn_mo_00_' res '.mat'];
    GRDname = ['MSK_' res '.mat'];
    % Tname = [dir0 'T_iocn_mo_01_60km_convt.mat'];
    % GRDname = ['MSK_240km_dst.mat'];
  otherwise
    disp('not a correct resolution');
end

% load MPAS GRID info
disp(['Now loading ' GRDname ' to: GRD']);
load([dir0 GRDname]); 

% load MPAS-O transport matrix
for i=1:length(TTM)
  Tname = sprintf('%s%s%s%s',TTM{i},fnp,res,'.mat'); 
  disp(['Now loading ' Tname]);
  load([dir0 Tname]);
  % ind = find(isnan(T)); T(ind)=0;
end



% TRdiv = -( A + D + H );

% load Lump and Spray operators
load([rdir 'Mask_TTM/Lump_Spray2.mat'],'L'); 
load([rdir 'Mask_TTM/Lump_Spray.mat'], 'S'); 

% A*h0 = 0;
%A = L*A*S;
%D = L*D*S;
%H = L*H*S;


P.spy          = 60*60*24*365;
P.tau          = (1/365)*P.spy;           % air-sea ideal age restoring timesale [s]

[ny,nz] = size(GRD.MASK);
iocn       = GRD.iocn;
n_iocn     = length(iocn);
P.Rmask    = GRD.MASK;

P.Ras = P.Rmask/P.tau;             % make air-sea rate operator (const)
P.Ras(:,2:end) = 0;
P.Ras = d0(P.Ras(iocn) );

h  = GRD.DZT(iocn);
%TRdiv = -d0(1./h)*( A +  H +  D*d0(h) ); 
% TRdiv = -(D);
TRdiv = -( d0(1./h)*H +  D  );

LHS = TRdiv +  P.Ras;
RHS = ones(length(iocn),1);

% steady state ideal age equation: TRdiv*a = 1
disp(sprintf('Factorize LHS of the %s equation',model));
% X = LHS \ RHS;
LP = mfactor(LHS);
X = mfactor(LP,RHS);

% xx = mat2mpas(X(end,:),GRD);
XM = mat2mpas(X,GRD); % to mpas-o grid (ny,nz)

% plot the 3th level of MPAS-O ideal age
%plotaa(GRD,XM,3,'v4');