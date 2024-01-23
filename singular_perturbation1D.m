function [A,uex,f]=singular_perturbation1D(n,ep)
% spectral differentiation (pseudospectral) matrix for 
% 1D singular pertubation equation with homogeneous Dirichlet boundary conditions
% -ep*\Delta u+u_x=f, 
% Input: n is the number of interior collocation nodes,
%        ep: epsilon
% Output: A is the associated  spectral differentiation matrix,
%         uex is the exact solution of the equation at the collocation
%         points.
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
e=sqrt(ep);
method=2; % 1=Legendre  2=Cheb
[D, node]=PSDirv(n,method);
X=node(2:n);
DD=D*D;
A=-ep*DD(2:n,2:n)+D(2:n,2:n);
%the exact solution and f for second order
uex=(X+1).*(1-exp((X-1)/e));
f=1+(2*e+X-(X+1)/e).*exp((X-1)/e);
%u=A\f;
%error=norm(uex-u,Inf), 
%plot(X,u,'-o',X,uex,'-')