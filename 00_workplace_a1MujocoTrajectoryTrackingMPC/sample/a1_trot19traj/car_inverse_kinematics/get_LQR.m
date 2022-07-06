clc
clear all
close all

dt = 0.2;
px = 0.05;
py = 0;
theta = pi/2;
A = eye(2);
c = cos(theta); s = sin(theta);

mat = [c-(py/px)*s, s+(py/px)*c;
         -(1/px)*s,    (1/px)*c]
     
B = mat*dt;

Q = eye(2);
R = 0.01*eye(2);
[K,S,E] = dlqr(A,B,Q,R)

K 
E
