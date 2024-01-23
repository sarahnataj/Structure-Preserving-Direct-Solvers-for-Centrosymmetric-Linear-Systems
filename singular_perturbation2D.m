function [A,uex,f]=singular_perturbation2D(n,ep)
% spectral differentiation (pseudospectral) matrix for 
% singular pertubation equation with homogeneous Dirichlet boundary conditions
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
[X,Y]=meshgrid(node(2:n), node(2:n));
DD=D*D;
A=-ep*kron(speye(n-1),DD(2:n,2:n))...
    -ep*kron(DD(2:n,2:n),speye(n-1))+kron(speye(n-1),D(2:n,2:n));
%the exact solution and f for second order
uex=(X+1).*(1-exp((X-1)/e)).*cos(pi*Y/2);
uex=reshape(uex',(n-1)^2,1);
f=ep*((2*exp((X - 1)/e).*(cos((pi*Y)/2))/e)+...
    (exp((X - 1)/e).*cos((pi*Y)/2).*((X + 1))/e^2) -...
    (pi^2*cos((pi*Y)/2).*((exp((X - 1)/e) - 1)).*((X + 1))/4)) - ...
    cos((pi*Y)/2).*((exp((X - 1)/e) - 1)) - ...
    (exp((X - 1)/e).*(cos((pi*Y)/2)).*((X + 1))/e);
f=reshape(f',(n-1)^2,1);
