% solve for the steady state radiocarbon in the ocean
clear
% load transport_Redi_Jan2013;
% set up transport and grid file
res = '60km';
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
disp(sprintf('load Transport operator T in %s',Tname));
load(Tname);
% load MPAS GRID info
disp(sprintf('load MPAS-O GRID GRD in %s',GRDname));
load(GRDname);

M3d = GRD.MASK;
TRdiv = -T;

sec_per_year = 365*24*60*60;
N = sum(M3d(:)); % total number of wet grid boxes
tau = 10*GRD.DZT(1,1)/50; % air-sea equilibration timescale
msk = M3d;
msk(:,:,2:end) = 0;
iocn = find(M3d(:)==1); % indeces of the wet points
L = d0(msk(iocn))/tau;
I = speye(N);
lam = log(2)/5730; % radiocarbon decay rate
T = TRdiv*sec_per_year;
Ro = ones(N,1);
b = L*Ro;
A = [T+L+lam*I];
% % LU decomposition of A
% tic
% FA = mfactor(A);
% toc
% tic
% R= mfactor(FA,b);
% toc
% %R = A\b;
% ILU and gmres to solve A\b
[X,fl0,rr0,it0,rv0] = myilu(A,b);
c14age = M3d+nan;
c14age(iocn) = -log(X)/lam;
save('c14age.mat','X','c14age','-v7.3');
