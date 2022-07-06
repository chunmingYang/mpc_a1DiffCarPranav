function zdot = diff_drive_rhs(t,z,u)
v = u(1); omega = u(2);
theta = z(3);
x_c_dot = v*cos(theta);
y_c_dot = v*sin(theta);
theta_dot = omega;
zdot = [x_c_dot y_c_dot theta_dot]';
