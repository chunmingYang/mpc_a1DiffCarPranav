clc
clear all
close all

t0 = 0;
tend = 480;
t = t0:0.01:tend;
tau = 2*pi*t/tend;
b = 2*pi/tend;
a = 3;

den = 1+sin(tau).*sin(tau);
 
x = a*cos(tau)./den; 
y = a*cos(tau).*sin(tau)./den;


xdot = (a*b*sin(tau).*(sin(tau).^2 - 3))./(sin(tau).^2 + 1).^2;
 
ydot = -(a*b*(3*sin(tau).^2 - 1))./(sin(tau).^2 + 1).^2;

figure(1)
subplot(2,1,1)
plot(t,x,'r'); hold on;
plot(t,y,'b');
legend('x','y');
subplot(2,1,2)
plot(t,xdot,'r'); hold on;
plot(t,ydot,'b');
legend('vx','vy');

% figure(2)
% plot(x,y);
% xlabel('x');
% ylabel('y');