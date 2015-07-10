%To simulate the system mx'' + kx = F
% z = [x;v]; {poisition; velocity]
close all;
param.m = 1;
param.k = 10;
param.omega = 0; % if omega is a non-zero value, it's a forced system

springmass = @(t,z)springprocess(z,t,param);

tspan = [0 10];
z0 = [0 1];

%used as a process model for the EKF
[t,z] = ode45(springmass,tspan,z0);

%Obtaining a noisy velocity measurement of known variance for the measurement model
variance = 1;
wn = sqrt(variance).*randn(size(z(:,2)));
y = z(:,2) + wn;


figure(1)
plot(t,z(:,1));
title('position vs time');
figure(2)
plot(t,z(:,2));
hold on
plot(t,y);
legend('state','measurement');
title('velocity vs time');
