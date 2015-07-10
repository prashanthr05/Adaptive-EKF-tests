%To simulate the system mx'' + kx = F
% z = [x;v]; {poisition; velocity]
close all;
param.m = 1;
param.k = 10;
param.omega = 2; % if omega is a non-zero value, it's a forced system

springmass = @(t,z)springprocess(z,t,param);

tspan = [0 10];
z0 = [0 1];

[t,z] = ode45(springmass,tspan,z0);
figure(1)
plot(t,z(:,1));
title('position vs time');
figure(2)
plot(t,z(:,2));
title('velocity vs time');
