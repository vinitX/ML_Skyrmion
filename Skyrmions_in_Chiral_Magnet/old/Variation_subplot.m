clear;
load('variation.mat');
f=zeros(size(flag,1),size(flag,2));
tol=0.1; % tolerance value for comparison
for i=1:size(fval,1)
for j=1:size(fval,2)
if flag(i,j)~=0 && fval(i,j)<0
    m1=vec_p(i,j,3);
    m2=vec_p(i,j,4);
    if abs(m1-m2)<tol && abs(m1)<tol && abs(m2)<tol
        f(i,j)=1;  ferromagnetic
    elseif m1>(1-tol) && m2>(1-tol)
        f(i,j)=2; % Spiral
    elseif abs(m1-m2)<tol && m1>-tol && m2>-tol
        f(i,j)=3; % Circular Cone
    elseif m1<m2 && m1>-tol && m2>-tol
        f(i,j)=4; % Elliptic Cone
    else
        f(i,j)=5; % Skyrmion or Others
    end
end
end
end
surf(f)
%surf(vec_p(:,:,1,1))
%surf(fval)