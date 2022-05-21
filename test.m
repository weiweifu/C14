% set up transport and grid file
rdir = '/DFS-L/DATA/moore/weiweif/MPAS-BGC/';
res = '240km';
% res = '60km';
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

% load MPAS-O transport matrix
for i=1:length(TTM)
  Tname = sprintf('%s%s%s%s',TTM{i},fnp,res,'.mat'); 
  disp(['Now loading ' Tname]);
  load([dir0 Tname]);
  % ind = find(isnan(T)); T(ind)=0;
end

% TRdiv = ( A + D + H );

% (H + A)*e1 = 0;
e = randn(size(A,1),1);
MF = mfactor(H + A);
for i = 1:20
  e1 = mfactor(MF,e);
  e1 = e1./(e1.'*e1);
end

h0 = e1;
TRdiv = ( A +  H + D.*spdiags(h0));
% steay state ideal age equation: TRdiv*a = -h0;
disp(sprintf('Factorize LHS of the %s equation',model));
% X = LHS \ RHS;
LP = mfactor(TRdiv);
X = mfactor(LP,-h0);

% load MPAS GRID info
disp(['Now loading ' GRDname ' to: GRD']);
load([dir0 GRDname]); 

% load Lump and Spray operators
load([rdir 'Mask_TTM/Lump_Spray2.mat'],'L'); 
load([rdir 'Mask_TTM/Lump_Spray.mat'], 'S'); 

P.spy          = 60*60*24*365;
P.lam          = log(2)/(5730*P.spy);   % c14 decay rate [s^-1]
if strcmp(model,'c14')
  P.tau          = 2*P.spy;               % air-sea c14 exchange timesale [s]
else
  P.tau          = 1/365*P.spy;           % air-sea ideal age restoring timesale [s]
end

[ny,nz] = size(GRD.MASK);
iocn       = GRD.iocn;
n_iocn     = length(iocn);
P.Rmask    = GRD.MASK;

P.R14 = P.lam*P.Rmask;             % make decay rate operator
P.R14 = d0(P.R14(iocn) );          

P.Ras = P.Rmask/P.tau;             % make air-sea rate operator (const)
P.Ras(:,2:end) = 0;
P.Ras = d0(P.Ras(iocn) );
P.So  = ones(n_iocn,1);            % atm C14 ratio (const)
 
P.init = P.Ras(iocn)*P.tau;

% $$ \frac{dR}{dt} + \left[\mathbf{T}+\lambda\mathbf{I}+\kappa\boldsymbol{\Lambda}\right]R = \kappa\boldsymbol{\Lambda}\mathbf{1}$$
if strcmp(model,'c14')
  disp(['Set up ' model ' model']);
  LHS = TRdiv + P.R14  + P.Ras;
  RHS = P.Ras*P.So;
else
  LHS = TRdiv +  P.Ras;
  RHS = ones(size(LHS,1),1);
end

% define function handle of RHS
if strcmp(model,'c14')
  fh = @(t,y) -LHS*y(:) + RHS;
else
  fh = @(t,y) -LHS*y(:) + RHS;
end

% disp(sprintf('Factorize LHS of the %s equation',model));
% X = LHS \ RHS;
% LP = mfactor(LHS);
% X = mfactor(LP,RHS)/P.spy;

% disp(sprintf('Time stepping (%s) : Euler backward',model));
% tspan = [ 0.0, P.spy*1.0/12/80 ];
% y0 = P.init;
% dt = 3600*0.5; 
% n = round(( tspan(2) - tspan(1) ) / dt);
% display(sprintf('There are %i iterations',n));
% [tE,X] = beuler_fixed ( fh, tspan, y0, n );

if strcmp(model,'c14')
  % c14age = -log(X)/(P.lam*P.spy);
  c14age = -log(X)/P.lam;
  XM = mat2mpas(c14age,GRD);
  save('c14age.mat','XM','X','c14age','-v7.3');
else
  % xx = mat2mpas(X(end,:),GRD);
  XM = mat2mpas(X,GRD); % to mpas-o grid (ny,nz)
  % save('iage.mat','X','XM','-v7.3');
end

% plot the 3th level of MPAS-O ideal age
plotaa(GRD,XM,3,'v4');
