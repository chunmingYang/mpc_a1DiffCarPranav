clc
clear all

syms t real
syms a b real

tau = b*t;


den = 1+sin(tau)*sin(tau);
x = a*cos(tau)/den; 
y = a*cos(tau)*sin(tau)/den;

xdot = simplify(diff(x,t))
ydot = simplify(diff(y,t))

xddot = simplify(diff(xdot,t))
yddot = simplify(diff(ydot,t))

%dydx = simplify(ydot/xdot)
%disp('dydx = ydot/xdot;');
%disp('theta = atan2(ydot,xdot);')