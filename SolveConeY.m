function y=SolveConeY(b,Y)
n=size(Y,1);
y=zeros(n,1);
m=floor(n/2);
if mod(n,2)==1
    y(m+1)=b(m+1)/Y(m+1,m+1);
end
for k=1:m
    p=m-k+1;
    q=m+k+mod(n,2);
    XX=[Y(p,p) Y(p,q);Y(p,q) Y(p,p)];
    if (p==m)&&(mod(n,2)==0)
        bb=[b(p)  b(q)]';
    else
        bb=[b(p)-Y(p,p+1:q-1)*y(p+1:q-1)  b(q)-Y(q,p+1:q-1)*y(p+1:q-1)]';
    end
    xx=XX\bb;
    y(p)=xx(1);y(q)=xx(2);
end
end