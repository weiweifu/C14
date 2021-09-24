function [x,varargout] = myilu(A,b,varargin)

disp('use ilu to generate the preconditioner');
[L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
tol = 1e-12;
maxit = 150;
disp(sprintf('use gmres to solve Ax=b: tol: %.0e maxit: %d',tol,maxit));
[x3,fl3,rr3,it3,rv3] = gmres(A,b,30,tol,maxit,L,U);
varargout{1} = fl3;
varargout{2} = rr3;
varargout{3} = it3;
varargout{4} = rv3;

end
