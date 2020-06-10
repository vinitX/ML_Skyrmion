%CFE is the function which optimizes the value of Continuum free energy by taking A, H, and an initial vector p0 as input.
function [p,fval] = CFE(A,H,p0)
%s is a function of p. three component of s corresponds to x, y and z axes respectively.
%p is set containing variables theta(i) and phi(i)
%the 21st variable of p contains L.
s=cell(3,10);
for i=1:10
    j=10+i;
    s{1,i}=@(p)(cos(p(j))*sin(p(i)));
    s{2,i}=@(p)(sin(p(j))*sin(p(i)));
    s{3,i}=@(p)cos(p(i));
end
F=@(p)0; %F is the array of Continum free energy functions

for i=1:9
    del=@(p)p(21)/10;%del is step size in L. del=L/100. 
    t1=@(p)(0.5*((s{1,i+1}(p)-s{1,i}(p))^2+(s{2,i+1}(p)-s{2,i}(p))^2+(s{3,i+1}(p)-s{3,i}(p))^2));% 0.5*grad(n)^2
    t2=@(p)(s{2,i}(p)*(s{3,i+1}(p)-s{3,i}(p))+s{3,i}(p)*(s{2,i+1}(p)-s{2,i}(p))); % n.curl(n)
    F=@(p)(F(p)+t1(p)*t1(p)/6-t1(p)+t2(p)+del(p)*(A*s{3,i}(p)^2-H*s{3,i}(p)));
end
[p,fval]=fminunc(F,p0); %Conjugate Gradiant Minimization
p(1:10)=mod(p(1:10),pi);
p(11:20)=mod(p(11:20),2*pi);