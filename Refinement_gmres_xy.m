% Given A, b in be given in precision u, this program solves 
% Ax = b using GMRES-based iterative refinement with inexact XY factorization 
% (calculated in single precision).( Algorithm 5).
%
% Input: A is the spectral differentiation (pseudospectral) matrix of size n
%       associated with
%            1D Poisson/variable--coefficient Poisson (diffusion) equations,
%            1D/2D biharmonic/variable--coefficient biharmonic equations,
%            3D Helmholtz equation.
%       u = the working precision at which A and b are stored=double,
%       uf = the precision at which the factorization of 𝐴 is computed = single,
%       ur = the precision at which the residual is calculated = double,
%       tol is the convergence criterion for the refinement process,
%       tolg is the tolerance for convergence of GMRES.
% Output: n_ir = the number of iterations for GMRES to converge,
%       e_u = norm of the error.
%
% Algorithm 5 is discussed in the manuscript: Structure-preserving 
% solvers for centrosymmetric linear systems with applications to spectral methods
% by Chen Greif, Sarah Nataj and Manfred Trummer.
%
% Author: Sarah Nataj, email:sarah.nataj@gmail.com
function [e_d, e_s, e_r,i, condition,condition_scaled,siz]=Refinement_gmres_xy(n,ecase)
%clc 
tol=1e-12;
tolg=1e-4;
scale=1;% use scaling to reduce the condition number
imax=100;
if ecase==1
[A,x_ext]=PseudoSpectral1D(n,4);%1D bihaarmonic
elseif ecase==2
[A,x_ext]=PseudoSpectral2D(n,4);%2D biharmonic
else
[A, C, M, x_ext]=Biharmonic_varcoef(n,1000);A=full(A);%2D variable coefficient biharmonic
%n=15;a=600;[A,x_ext]=Helmholtz3D(n,a);A=full(A);
end
siz=size(A,1);
condition=cond(A);
b=A*x_ext;
maxit = size(A,1);
if scale==1
    [R,S]=scaling(A);
    Ah=R*A*S;bh=R*b;
    condition_scaled=cond(Ah);
else
    Ah=A;bh=b;
end

%double
[Pd,Xd,Yd]=xy(Ah,0);
w=SolveConeXd(Xd,Pd*bh);
y=SolveConeYd(Yd,w);
if scale==1
x_u=S*y;
else
 x_u=y;
end
e_d = (x_ext-x_u)/x_ext;
e_d = norm(e_d);


%single
[P,X,Y]=xy(single(Ah),0);
w=SolveConeXs(X,P*bh);
y=SolveConeYs(Y,w);
if scale==1
x_u=double(S*y);
else
 x_u=double(y);
end
e_s = double(x_ext-x_u)/x_ext;
e_s = double(norm(e_s));

% refinement with single precision factors in GMRES
% solve the residual equation
for i=1:imax
r_u = double(b-A*x_u);
if scale==1
    r_u=R*r_u;
    [z_u,fl,rr,it,rv]=gmres(P*Ah,P*r_u,[],tolg,maxit,X,Y);z_u=S*z_u;
else
    [z_u,fl,rr,it,rv]=gmres(P*A,P*r_u,[],tolg,maxit,X,Y);
end
x_u=x_u+z_u;
e_r = double(x_ext-x_u)/x_ext;
e_r = double(norm(e_r));
if e_r<=tol
    %x_u;e_u,i
    return 
end

end

end

function x=SolveConeXs(X,b)
X=single(X);b=single(b);
n=size(X,1);
x=single(zeros(n,1));
m=floor(n/2);
for k=1:m
    p=k;
    q=n-k+1;
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==1)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,1:p-1)*x(1:p-1)-X(p,q+1:end)*x(q+1:end),...
            b(q)-X(q,1:p-1)*x(1:p-1)-X(q,q+1:end)*x(q+1:end)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
if mod(n,2)==1
    x(m+1)=single((b(m+1)-X(m+1,:)*x)/X(m+1,m+1));
end
end

function x=SolveConeYs(X,b)
X=single(X);b=single(b);
n=size(X,1);
x=single(zeros(n,1));
m=floor(n/2);
if mod(n,2)==1
    x(m+1)=b(m+1)/X(m+1,m+1);
end
for k=1:m
    p=m-k+1;
    q=m+k+mod(n,2);
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==m)&&(mod(n,2)==0)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,p+1:q-1)*x(p+1:q-1)  b(q)-X(q,p+1:q-1)*x(p+1:q-1)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
end

function x=SolveConeXd(X,b)
n=size(X,1);
x=zeros(n,1);
m=floor(n/2);
for k=1:m
    p=k;
    q=n-k+1;
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==1)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,1:p-1)*x(1:p-1)-X(p,q+1:end)*x(q+1:end),...
            b(q)-X(q,1:p-1)*x(1:p-1)-X(q,q+1:end)*x(q+1:end)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
if mod(n,2)==1
    x(m+1)=(b(m+1)-X(m+1,:)*x)/X(m+1,m+1);
end
end

function x=SolveConeYd(X,b)
n=size(X,1);
x=zeros(n,1);
m=floor(n/2);
if mod(n,2)==1
    x(m+1)=b(m+1)/X(m+1,m+1);
end
for k=1:m
    p=m-k+1;
    q=m+k+mod(n,2);
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==m)&&(mod(n,2)==0)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,p+1:q-1)*x(p+1:q-1)  b(q)-X(q,p+1:q-1)*x(p+1:q-1)]';
    end
    xx=XX\bb;
    x(p)=xx(1);x(q)=xx(2);
end
end



