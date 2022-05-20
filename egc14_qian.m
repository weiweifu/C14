% solve for the steady state radiocarbon in the ocean
clear
load transport_Redi_Jan2013;
sec_per_year = 365*24*60*60;
N = sum(M3d(:)); % total number of wet grid boxes
tau = 10*grd.dzt(1)/50; % air-sea equilibration timescale
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
% LU decomposition of A
tic
FA = mfactor(A);
toc
tic
R= mfactor(FA,b);
toc
%R = A\b;
c14age = M3d+nan;
c14age(iocn) = -log(R)/lam;
