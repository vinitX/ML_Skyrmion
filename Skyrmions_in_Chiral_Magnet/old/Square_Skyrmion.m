clear;
K=1.5;
H=0;
D=1.5;
%s is a function of p. three component of s corresponds to x, y and z axes respectively.
%p is set containing variables theta(i,j) and phi(i,j+size)
size=4;
L=size;
del=L/size;
p0=[size,2*size];
for i=1:size
    for j=1:size
        r=(i^2+j^2)^0.5;
        p0(i,j)=pi*r/L;
        p0(i,j+size)=asin(i/r);
    end
end

s=cell(3,size+1,size+1);
for i=1:size+1
    for j=1:size+1
        l=i;
        m=j;
        if i>size  %boundary conditions
                l=l-size;
        end
        if j>size
                m=m-size;
        end
        k=size+m;
        s{1,i,j}=@(p)(cos(p(l,k))*sin(p(l,m)))*del;
        s{2,i,j}=@(p)(sin(p(l,k))*sin(p(l,m)))*del;
        s{3,i,j}=@(p)cos(p(l,m))*del;
    end
end

F=@(p)0; %F is the array of Continum free energy functions
for i=1:size
    for j=1:size
        v1=@(p)0.5*((s{1,i+1,j}(p)-s{1,i,j}(p))^2+(s{2,i+1,j}(p)-s{2,i,j}(p))^2+(s{3,i+1,j}(p)-s{3,i,j}(p))^2)/(del^2);
        v2=@(p)0.5*((s{1,i,j+1}(p)-s{1,i,j}(p))^2+(s{2,i,j+1}(p)-s{2,i,j}(p))^2+(s{3,i,j+1}(p)-s{3,i,j}(p))^2)/(del^2);
        t1=@(p)v1(p)+v2(p); %0.5*grad(n)^2
        u1=@(p)(s{3,i,j+1}(p)-s{3,i,j}(p))/del;
        u2=@(p)(s{3,i+1,j}(p)-s{3,i,j}(p))/del;
        u3=@(p)(s{2,i+1,j}(p)-s{2,i,j}(p))/del;
        u4=@(p)(s{1,i,j+1}(p)-s{1,i,j}(p))/del;
        t2=@(p)(s{1,i,j}(p)*(u1(p))-s{2,i,j}(p)*(u2(p))+s{3,i,j}(p)*(u3(p)-u4(p))); % n.curl(n)
        F=@(p)del*(F(p)+t1(p)^2/6-t1(p)+D*t2(p)-(K/2)*(1-s{3,i,j}(p)^2)+(H/4)*(1-s{3,i,j}(p)));
    end
end
[p,fval]=fminunc(F,p0); %Conjugate Gradiant Minimization
save('square.mat','p','fval');

theta=mod(p(:,1:size),pi);
phi=mod(p(:,size+1:2*size),2*pi);
%x, y and z are the components of n(x,y) along i, j and k directions.
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta);

%R=[cos(b) sin(b) 0 ; -sin(b) cos(b) 0 ; 0 0 1];
X=ones(size,size);
quiver3(X,x,y,z);