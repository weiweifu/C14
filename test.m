% set up transport and grid file
res = '240km';
% mach = 'cori';
mach = 'gp';

switch lower(mach)
  case 'cori'
    dir0 = '/global/cfs/cdirs/e3sm/weiweif/MPAS_IRF/Mask_TTM/';
    Tname = [dir0 'T_iocn_mo_00_' res '.mat'];
    GRDname = ['MSK_' res '.mat'];
  case 'gp'
    dir0 = '/DFS-L/DATA/moore/weiweif/MPAS-BGC/Data/';
    Tname = [dir0 'T_iocn_mo_00_' res '.mat'];
    GRDname = [dir0 'MSK_' res '.mat'];
  otherwise
    disp('not a correct resolution');
end

% load MPAS-O transport matrix
disp('load Transport operator: T');
load(Tname); 
% load MPAS GRID info
disp('load MPAS-O GRID: GRD');
load(GRDname); 

P.spyr         = 60*60*24*365;
P.lam          = log(2)/(5730*P.spyr);   % c14 decay rate [s^-1]
P.tau          = 2*P.spyr;               % air-sea c14 exchange timesale [s]

[ny,nz] = size(GRD.MASK);
iocn       = GRD.iocn;
n_iocn     = length(iocn);
P.Rmask    = GRD.MASK;

P.R14 = P.lam*P.Rmask;             % make decay rate operator
P.R14 = d0(P.R14(iocn) );          

P.Ras = P.Rmask/P.tau;             % make air-sea rate operator (const)
P.Ras(:,:,2:end) = 0;
P.Ras = d0(P.Ras(iocn) );
P.So  = ones(n_iocn,1);            % atm C14 ratio (const)
 
% $$ \frac{dR}{dt} + \left[\mathbf{T}+\lambda\mathbf{I}+\kappa\boldsymbol{\Lambda}\right]R = \kappa\boldsymbol{\Lambda}\mathbf{1}$$
LHS = T + P.R14  + P.Ras;
RHS = P.Ras*P.So;

disp('Factorize LHS of the C14 equation');
% X = LHS \ RHS;
LP = mfactor(LHS);
X = mfactor(LP,RHS);
c14age = -log(X)/(P.lam*P.spyr);
save('c14age.mat','X','c14age','-v7.3');
