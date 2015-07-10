function dzdt = springprocess(z,t,p)

k = p.k;
m = p.m;
omega = p.omega;


dx = z(2);



dv = (-k*z(1) + 10*sin(omega*t))/m;

% dv = (-k*z(1))/m;
dzdt = [dx;dv];
end