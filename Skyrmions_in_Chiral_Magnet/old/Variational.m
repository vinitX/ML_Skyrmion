clear;
%p is a vector with seven variables corresponding to q, m and t in order.
q0=[1 1];
m0=[0.5 0.5];
t0=[pi/2 pi/2 pi/2];
p0=[q0 m0 t0]; 
b=pi/2; %beta
J1=1;
J2=1;
D=0:0.1:3;
K=-2:0.1:4;
H=-2:0.1:4;

fval=zeros(length(D),length(K),length(H));
vec_p=zeros(length(D),length(K),length(H),7);
flag=zeros(length(D),length(K),length(H));

for i=1:length(D)
    for j=1:length(K)
        for k=1:length(H)
        fj1=@(p)(-J1/8)*(p(1)^2+p(2)^2)*(1-((1-p(3)^2)*(1-p(4)^2))^0.5);
        fj2=@(p)(J2/96)*(p(1)^2+p(2)^2)*(4-2*(2+p(3)^2+p(4)^2)*((1-p(3)^2)*(1-p(4)^2))^0.5);
        fd=@(p)-D(i)*p(1)*p(3)*p(4)*sin(p(5))*sin(p(6));
        fk=@(p)(K(j)/2)*cos(p(6))^2*(2-p(3)^2-p(4)^2)+(K(j)/2)*sin(p(6))^2*((p(3)*cos(p(7)))^2+(p(4)*sin(p(7)))^2);
        %This function computes Elliptic Integral of Second Kind.
        step=(pi/2)/10;
        ellip2=@(p)0;
        for theta=0:step:pi/2
            ellip2=@(p)ellip2(p)+(1-p*(sin(theta)^2))*step;
        end
        E=@(p)ellip2(-(p(3)^2-p(4)^2)/(1-p(3)^2));
        fh=@(p)(-2*H(k)/pi)*((1-p(3)^2)^0.5)*cos(p(6))*E(p);
        F=@(p)fj1(p)+fj2(p)+fd(p)+fk(p)+fh(p);
        [vec_p(i,j,k,:),fval(i,j,k),flag(i,j,k)]=fminunc(F,p0);
        end
    end
end
save('variation.mat','vec_p','fval','flag');