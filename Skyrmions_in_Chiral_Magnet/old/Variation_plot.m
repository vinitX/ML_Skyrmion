clear;
load('variation (2).mat');
tol=0.15; % tolerance value for comparison
for i=1:size(fval,1)
for j=1:size(fval,2)
for k=1:size(fval,3)
if flag(i,j,k)~=0
if abs(fval(i,j,k))<10^4 && fval(i,j,k)<0
    m1=vec_p(i,j,k,3);
    m2=vec_p(i,j,k,4);
    if abs(m1-m2)<tol && abs(m1)<tol && abs(m2)<tol
        flag(i,j,k)=1; %#ok<*SAGROW> % ferromagnetic
    elseif m1>(1-tol) && m2>(1-tol)
        flag(i,j,k)=2; % Spiral
    elseif abs(m1-m2)<tol && m1>-tol && m2>-tol
        flag(i,j,k)=3; % Circular Cone
    elseif m1<m2 && m1>-tol && m2>-tol
        flag(i,j,k)=4; % Elliptic Cone
    else
        flag(i,j,k)=5; % Skyrmion
    end
else
    flag(i,j,k)=-1; % Skyrmion or Nothing
end
end
end
end
end
surf(reshape(flag(:,:,40),11,61))