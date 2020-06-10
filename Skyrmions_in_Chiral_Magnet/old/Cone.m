%Cone_plot takes the value of A and H and which optimizes the value of 
%Continuum free energy using a initial vector p0. and plots the data to
%visualize a elliptic cone in Dresselhaus Limit.
clear;
K=0.5;
H=0;
D=0;

%s is a function of p. three component of s corresponds to x, y and z axes respectively.
%p is set containing variables theta(i) and phi(i)
%the 21st variable of p contains L.
size=8;
p0=ones(1,2*size+1);
p0(2*size+1)=size;
for i=1:size
    p0(i)=pi/2;
    p0(i+size)=i*2*pi/size;
end
del=@(p)abs(p(2*size+1)/size);%del is step size in L. del=L/size.

s=cell(3,size);
for i=1:size
    j=size+i;
    s{1,i}=@(p)(cos(p(j))*sin(p(i)))*del(p);
    s{2,i}=@(p)(sin(p(j))*sin(p(i)))*del(p);
    s{3,i}=@(p)cos(p(i))*del(p);
end
%Boundary Conditions
s{1,size+1}=@(p)(cos(p(size+1))*sin(p(1)))*del(p);
s{2,size+1}=@(p)(sin(p(size+1))*sin(p(1)))*del(p);
s{3,size+1}=@(p)cos(p(1))*del(p);

F=@(p)0; %F is the Continum free energy functions

for i=1:size
    t1=@(p)0.5*((s{1,i+1}(p)-s{1,i}(p))^2+(s{2,i+1}(p)-s{2,i}(p))^2+(s{3,i+1}(p)-s{3,i}(p))^2)/(del(p)^2);% 0.5*grad(n)^2
    t2=@(p)(s{2,i}(p)*(s{1,i+1}(p)-s{1,i}(p))+s{1,i}(p)*(s{2,i+1}(p)-s{2,i}(p)))/del(p); % n.curl(n)
    F=@(p)F(p)+del(p)*(t1(p)^2-t1(p)+D*t2(p)-(K/2)*(1-s{3,i}(p)^2)+(H/4)*(1-s{3,i}(p)));
end
[p,fval]=fminunc(F,p0); %Conjugate Gradient Minimization.
theta=mod(p(1:size),pi);
phi=mod(p(size+1:2*size),2*pi);
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta);
X=zeros(1,size);
Y=zeros(1,size);
Z=1:size;
quiver3(X,Y,Z,x,y,z);
