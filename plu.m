function [L, U, P]=plu(A)
n=size(A,1);
p=1:n;
for k=1:n-1
    [val,q]=max(abs(A(k:n,k)));
    q=q+(k-1);
    A([k,q],:)=A([q,k],:);
    p([k,q])=p([q,k]);
    J=k+1:n;
    A(J,k)=A(J,k)/A(k,k);
    A(J,J)=A(J,J)-A(J,k)*A(k,J);
end
L=zeros(n,n);U=zeros(n,n);
for k=1:n
    L(k,1:k-1)=A(k,1:k-1);L(k,k)=1;
    U(k,k:n)=A(k,k:n);
    I=eye(n);
    P=I(p,:);
end
end
