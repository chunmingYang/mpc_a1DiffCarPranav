function diff_drive_ik
clc
clear 
close all

%%%%%% parm is struct that helps pass many argument %%
parms.R = 0.2; %length of the robot
parms.px = 0.05; %0.05;  %x-location of the pen from mid-point
parms.py = 0;        %y-location of the pen from mid-point
parms.Kp = 1; %gain for controller
a = 3;

parms.writeMovie = 0; %set to 1 to get a movie output
parms.nameMovie = 'car_drawing.avi';

t0 = 0;
tend = 120;
fps = 2; %adjust number of frames per second
parms.delay = 0.001; %adjust the delay in the animation
t = t0:0.01:tend;

%%% circle %%%
%x_ref = 1*cos(2*pi*t/tend); y_ref = 1*sin(2*pi*t/tend);

% %%% astroid (note these are coordinates of point P) %%
% x_center = 0; y_center = 0; a = 5;
% x_ref = x_center+a*cos(2*pi*t/tend).^3; y_ref = y_center+a*sin(2*pi*t/tend).^3;

% %% leminscate
x_center = -a; y_center = 0;
%tau = 2*pi*t/tend;
b = 2*pi/tend;
tau = b*t;
den = 1+sin(tau).*sin(tau);
den2 = den.*den;
x_ref = x_center+a*cos(tau)./den; y_ref = y_center+a*cos(tau).*sin(tau)./den;
xdot_ref = (a*b*sin(tau).*(sin(tau).*sin(tau) - 3))./den2;
ydot_ref = -(a*b*(3*sin(tau).*sin(tau) - 1))./den2;

theta0 = pi/2; %guessed to be zero. Can be made better by specifying a better value based on the specific curve
[x0,y0] = ptP_to_ptC(x_ref(1),y_ref(1),theta0,parms); %given the x_ref(1) and y_ref(1) (position of E at start) this step computes the correct starting position for the car (i.e., point P)
z0 = [x0 y0 theta0]; %x, y, theta at the start for point C
z = z0;
e = [0 0];
v = 0;
omega = 0;

for i=1:length(t)-1
    
    % 1. get x_c,y_x position
    x_c = z(end,1); y_c = z(end,2); theta = z(end,3);
    
    % 2. get x_p, y_p from x_c,y_c 
    [x_p,y_p] = ptC_to_ptP(x_c,y_c,theta,parms);
    
    % 3. get error = xe-x_ref and ye-y_ref
    error = [x_ref(i)-x_p y_ref(i)-y_p];
    e = [e; error];
   
    %4. get u = [v, omega] from the errors
     err_p_dot = [xdot_ref(i); ydot_ref(i)] + [parms.Kp*error(1); parms.Kp*error(2)];
     c = cos(theta); s = sin(theta);
     px = parms.px; py = parms.py;
     A = [c-(py/px)*s s+(py/px)*c; ...
           -(1/px)*s     (1/px)*c];
     u = A*err_p_dot; %u = [v omega]
     v = [v; u(1)]; %save v
     omega = [omega; u(2)]; %save omega

    
    % 5. now control the car based on u = [v omega]
    % We use a matlab integrator ode4 which is based on runge-kutta method
    zz = ode4(@diff_drive_rhs,[t(i) t(i+1)],z0,u);
    z0 = zz(end,:);
    z = [z; z0];
end


%%%%%% get coarse data for animation
t_interp = linspace(t0,tend,fps*tend);
[m,n] = size(z);
for i=1:n
    z_interp(:,i) = interp1(t,z(:,i),t_interp);
end

% figure(1)
% animation(t_interp,z_interp,parms);


figure(2)
subplot(2,1,1)
plot(t,v,'m'); hold on
plot(t,omega,'c');
ylabel('Velocity');
legend('v','\omega','Location','Best');
subplot(2,1,2)
plot(t,e(:,1),'r'); hold on
plot(t,e(:,2),'b'); 
legend('error x','error y','Location','Best');
ylabel('error');
xlabel('time');


function [x_p,y_p] = ptC_to_ptP(x_c,y_c,theta,parms)
c = cos(theta); s = sin(theta); 
p = [parms.px; parms.py];
Xp = [c -s; s c]*p + [x_c; y_c];
x_p = Xp(1); y_p = Xp(2);

function [x_c,y_c] = ptP_to_ptC(x_p,y_p,theta,parms)
c = cos(theta); s = sin(theta); 
p = [parms.px; parms.py];
Xc = -[c -s; s c]*p + [x_p; y_p];
x_c = Xc(1); y_c = Xc(2);