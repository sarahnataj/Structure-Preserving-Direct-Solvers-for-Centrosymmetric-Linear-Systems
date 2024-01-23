function y=SolveConeX(b,X)
n=size(X,1);
y=zeros(n,1);
m=floor(n/2);
for k=1:m
    p=k;
    q=n-k+1;
    XX=[X(p,p) X(p,q);X(p,q) X(p,p)];
    if (p==1)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-X(p,1:p-1)*y(1:p-1)-X(p,q+1:end)*y(q+1:end)  b(q)-X(q,1:p-1)*y(1:p-1)-X(q,q+1:end)*y(q+1:end)]';
    end
    xx=XX\bb;
    y(p)=xx(1);y(q)=xx(2);
end
if mod(n,2)==1
    y(m+1)=(b(m+1)-X(m+1,:)*y)/X(m+1,m+1);
end
end