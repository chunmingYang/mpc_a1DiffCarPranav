
clc
clear all
close all

epsilon = 1e-6;
a = 2;
t0 = 0; tf = 90;
b = 2*pi/(tf-t0);

%t = linspace(t0,tf);
t = 0:0.2:90;
tau = b*t;

den = 1+sin(tau).*sin(tau);
x = a*cos(tau)./den; 
y = a*cos(tau).*sin(tau)./den;


xdot = (a*b*sin(tau).*(sin(b*t).^2 - 3))./(sin(tau).^2 + 1).^2;
ydot = (a*b*(3*cos(2*tau) - 1))./(2*(sin(tau).^2 + 1).^2);

 
xddot = (a*b^2*cos(tau).*(10*cos(tau).^2 + cos(tau).^4 - 8))./(cos(tau).^2 - 2).^3;
yddot = -(2*a*b^2*cos(tau).*sin(tau).*(3*cos(tau).^2 + 2))./(sin(tau).^2 + 1).^3;
 
theta = atan2(ydot,xdot);
sec_theta = sec(theta);
num2 = xdot.*yddot - ydot.*xddot;
den2 = (xdot.*xdot+epsilon);
thetadot = (1./sec_theta.^2).*(num2./den2);

theta_int = zeros(1,length(theta));
theta_int(1) = theta(1);
%thetadot(1) = 0;
for i=1:length(theta_int)-1
    theta_int(i+1) = theta_int(i) + (t(i+1)-t(i))*thetadot(i);
end
%theta_int


fig_no  = 1;

figure(fig_no); fig_no = fig_no + 1;
subplot(2,1,1);
plot(t,x,'r'); hold on;
plot(t,y,'b-.');
subplot(2,1,2);
plot(t,xdot,'r'); hold on;
plot(t,ydot,'b-.');

figure(fig_no); fig_no = fig_no + 1;
plot(x,y,'r','Linewidth',2);

figure(fig_no); fig_no = fig_no + 1;
subplot(2,1,1)
plot(t,theta,'k.','Linewidth',2); hold on
plot(t,theta_int,'r','Linewidth',2);
% plot(t,-cos(theta),'b','Linewidth',2);
subplot(2,1,2);
plot(t,thetadot,'b','Linewidth',2);

