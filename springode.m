%To simulate the system mx'' + kx = F
% z = [x;v]; {position; velocity]
close all;
clear;
clc

%% Adding paths
utilities = genpath('./utils');
addpath(utilities);

%% Initialize process variables
param.m = 1;
param.k = 1;
param.omega = 0; % if omega is a non-zero value, it's a forced system
param.Amp = 0; %Amplitude for forced input

Ts = 0.001;
% Initial conditions for the ODE solver
tspan = 0:Ts:10; %[0 10];
z0 = [2 0];
%% Use ODE solver to integrate the system through time
springmass = @(t,z)springprocess(z,t,param);
%used as a reference for process model for the EKF
[t,z] = ode45(springmass,tspan,z0);

%% Obtaining a noisy velocity measurement of known variance for the measurement model
variance = 0.50;
wn = sqrt(variance).*randn(size(z(:,2))); % white noise with a given variance

y = z(:,2) + wn;

%% Initializing the filter params
n = 2; %state dimension
m = 1; %output dimension
%% Define the coefficient matrices
A_cont= [0 1;...
    -param.k/param.m 0];      %State Transition matrix
B_cont = [0; 1/param.m];           % Input effect matrix
C_cont = [0 1]; 

%Ts= 0.001
A = expm(A_cont*Ts);
B = inv(A_cont)*(A - eye(2))*B_cont;
C = C_cont;

%% Initialize the covariance matrices
P = 0.01*ones(n,n);%2*eye(n); % NxN state covariance matrix
Q = zeros(n,n);%25*eye(n);    % randomly assumed Process noise covariance
R = variance;
 
%% Initialize the filter state estimate and covariance
xh = z0'; %Initial x0
Ph = P;


%% Run the Filter
for i = 1:size(t,1)
    U = param.Amp*sin(param.omega*t(i)); %sinusoidal input
%     U = param.Amp*t(i)*t(i);  %ramp input
%     U = param.Amp*t(i);  %step input

% KF update
    [xupdt Pupdt] = kf_update(xh,Ph,y(i),C,R);
    xUpt(i,:) = xupdt;
    xh = xupdt;
    Ph = Pupdt;
    Ph = (Ph + Ph')/2;

    
% KF predict
    [xhat Phat] =  kf_predict(xh,Ph,A,Q,B,U);
    xPred(i,:) = xhat;
    xh = xhat;
    Ph = Phat;
    %
end
%
%
%% Plot the desired behavior, measurement, predicted and corrected estimates
figure(1)
plot(t,z(:,1),'b');
hold on
pause;
plot(t,xPred(:,1),'r');
hold on
pause;
plot(t,xUpt(:,1),'g');
legend('actual','predicted','corrected');
title('position vs time');
% legend('actual','measurement','predicted','corrected');

figure(2)
plot(t,z(:,2),'g');
title('velocity vs time');
hold on
pause;
plot(t,y,'y');
hold on
pause;
plot(t,xPred(:,2),'r');
hold on
pause;
plot(t,xUpt(:,2),'k');
legend('actual','measurement','predicted','corrected');
% legend('actual','predicted','corrected');
