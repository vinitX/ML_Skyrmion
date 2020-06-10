clear;
H=0:1:4.1;
A=-1.6:1:3.1;
p0=ones(1,21); %intial vector p0. theta=pi/2. 
p0(1:10)=pi/2;
p0(11:20)=rand(1,10)*2*pi;

F=zeros(5,5); %set of optimized free energy for different values of H and A.
p=zeros(5,5,21); %p stores the optimized values of theta, phi and L for different values of H and A.
for i=1:5
    for j=1:5
        [p(i,j,:),F(i,j)]=CFE(A(i),H(j),p0);
    end
end
