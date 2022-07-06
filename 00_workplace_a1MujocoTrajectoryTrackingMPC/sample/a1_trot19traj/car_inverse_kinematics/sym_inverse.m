syms cx cy theta real
% A = [cos(theta) -px*sin(theta)-py*cos(theta); ...
%      sin(theta)  px*cos(theta)-py*sin(theta)];
%  
%  simplify(inv(A))
 
c = [cx; cy];
A = [rotation(theta) rotation(pi/2+theta)*c; 0 0 1];
simplify(A)
simplify(inv(A))

function R = rotation(theta)

c = cos(theta); s = sin(theta);
R = [c -s; s c];
end