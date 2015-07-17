function dzdt = springprocess(z,t,p)

k = p.k;
m = p.m;
omega = p.omega;
Amp = p.Amp;

dx = z(2);
% dv = (-k*z(1)/m); %no input
dv = (-k*z(1)-sign(z(1))*k*z(1)*z(1) + Amp*sin(omega*t))/m; %sinusoidal input
% dv = (-k*z(1) + Amp*t*t)/m;          %ramp input
% dv = (-k*z(1) + Amp*t)/m;            %step input


dzdt = [dx;dv];
end