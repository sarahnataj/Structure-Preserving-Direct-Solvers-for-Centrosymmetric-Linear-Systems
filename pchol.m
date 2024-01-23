function L=pchol(A)
n=size(A,1);
for k=1:n-1
    A(k,k)=sqrt(A(k,k));
    for i=k+1:n
        A(i,k)=A(i,k)/A(k,k);
    end
    for j=k+1:n
        for i=j:n
            A(i,j)=A(i,j)-A(i,k)*A(j,k);
        end
    end
end
A(n,n)=sqrt(A(n,n));
L=zeros(n,n);
for k=1:n
    L(k,1:k)=A(k,1:k);
end
L=L';
end