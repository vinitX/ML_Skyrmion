clear;
%p is a vector with seven variables corresponding to q, m and t in order.
q0=[0 0];
m0=[0 0];
t0=[0 0 0];
p0=[q0 m0 t0]; 
b=pi/2; %beta
J1=1;
J2=1;
D=0:0.05:4;
K=0:0.05:4;
H=0;

fval=zeros(length(D),length(K),length(H));
vec_p=zeros(length(D),length(K),length(H),7);
flag=zeros(length(D),length(K),length(H));

step=(pi/2)/1000;
ellip2=@(p)0;
for theta=0:step:pi/2
    ellip2=@(p)ellip2(p)+((1-p*(sin(theta)^2))^0.5)*step;
end
        
for i=1:length(D)
    for j=1:length(K)
        for k=1:length(H)
        fj1=@(p)(-J1/8)*(p(1)^2+p(2)^2)*(1-((1-p(3)^2)*(1-p(4)^2))^0.5);
        fj2=@(p)(J2/96)*(p(1)^2+p(2)^2)*(4-2*(2+p(3)^2+p(4)^2)*((1-p(3)^2)*(1-p(4)^2))^0.5);
        fd=@(p)-D(i)*p(1)*p(3)*p(4)*sin(p(5))*sin(p(6));
        fk=@(p)(K(j)/2)*cos(p(6))^2*(2-p(3)^2-p(4)^2)+(K(j)/2)*sin(p(6))^2*((p(3)*cos(p(7)))^2+(p(4)*sin(p(7)))^2);
        %This function computes Elliptic Integral of Second Kind.
        
        E=@(p)ellip2(-(p(3)^2-p(4)^2)/(1-p(3)^2));
        fh=@(p)(-2*H(k)/pi)*((1-p(3)^2)^0.5)*cos(p(6))*E(p);
        F=@(p)fj1(p)+fj2(p)+fd(p)+fk(p)+fh(p);
        lb=[-1,-1,-1,-1,-pi,-pi,-pi];
        ub=[1,1,1,1,pi,pi,pi];
        [vec_p(i,j,k,:),fval(i,j,k),flag(i,j,k)]=fmincon(F,p0,[],[],[],[],lb,ub);
        end
    end
end
save('variation.mat','vec_p','fval','flag');