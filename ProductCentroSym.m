function C=ProductCentroSym(A,B)
%Find product of two centrosymmetric matrices, A is m*p and B is p*n centrosymmetric matrices. 
%A = centro_generator(6, 6);
%B = centro_generator(6, 6);
m=size(A,1);
n=size(B,2);
M=A(1:floor(m/2),:)*B;
M=reshape(M',1,[]);
if mod(m,2) == 1
    c=A(ceil(m/2),:)*B(:,1:ceil(n/2));
    M = [M, c];
    if mod(n,2)==1
    M = [M, fliplr(M(1:(end-1)))];
    else
        M = [M, fliplr(M)];
    end
    C = reshape(M, [n m])';
else
    M = [M, fliplr(M)];
    C= reshape(M, [n m])';
end
%C,A*B, norm(C-A*B)
end