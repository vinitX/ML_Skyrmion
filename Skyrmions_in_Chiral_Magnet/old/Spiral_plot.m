%Spiral_plot takes the value of A and H and which optimizes the value of 
%Continuum free energy using a initial vector p0. and plots the data to
%visualize a spiral in Dresselhaus Limit.
clear;
A=0;
H=0;
%s is a function of p. three component of s corresponds to x, y and z axes respectively.
%p is set containing variables theta(i) and phi(i)
%the 21st variable of p contains L.
size=15;
p0=[mod(rand(1,size),pi) mod(rand(1,size),2*pi) rand(1,1)];
s=cell(3,size);
for i=1:size
    j=size+i;
    s{1,i}=@(p)(cos(p(j))*sin(p(i)));
    s{2,i}=@(p)(sin(p(j))*sin(p(i)));
    s{3,i}=@(p)cos(p(i));
end
%Boundary Conditions
s{1,size+1}=@(p)(cos(p(size+1))*sin(p(1)));
s{2,size+1}=@(p)(sin(p(size+1))*sin(p(1)));
s{3,size+1}=@(p)cos(p(1));

F=@(p)0; %F is the array of Continum free energy functions

for i=1:size-1
    del=@(p)p(21)/10;%del is step size in L. del=L/100. 
    t1=@(p)(0.5*((s{1,i+1}(p)-s{1,i}(p))^2+(s{2,i+1}(p)-s{2,i}(p))^2+(s{3,i+1}(p)-s{3,i}(p))^2))/(del(p)^2);% 0.5*grad(n)^2
    t2=@(p)(s{2,i}(p)*(s{3,i+1}(p)-s{3,i}(p))+s{3,i}(p)*(s{2,i+1}(p)-s{2,i}(p)))/del(p); % n.curl(n)
    F=@(p)(F(p)+t1(p)*t1(p)/6-t1(p)+t2(p)+(A*s{3,i}(p)^2-H*s{3,i}(p)));
end
[p,fval]=fminunc(F,p0); %Conjugate Gradiant Minimization
theta=mod(p(1:size),pi);
phi=mod(p(size+1:2*size),2*pi);
x=cos(phi).*sin(theta);
y=sin(phi).*sin(theta);
z=cos(theta);
X=1:size;
quiver3(X,x,y,z);