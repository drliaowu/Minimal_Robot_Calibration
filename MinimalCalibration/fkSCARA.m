function [ T ] = fkSCARA(xi,theta,n)
%[ T ] = fkSCARA(xi,theta,n) : forward kinematics of Scara type robot
%
%Input:
%   xi=6*n, joint twists;
%   theta=n*1, joint viariables;
%   n, number of joints
%
%Output:
%   T=4*4, homogeneous transformation

T=eye(4);
for i=1:n
    T=T*se3Exp(xi(:,i)*theta(i));
end
end

