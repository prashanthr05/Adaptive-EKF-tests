%To simulate the system mx'' + kx = F
% z = [x;v]; {position; velocity]
close all;
clear;
clc

%% Adding paths
utilities = genpath('./utils');
dynFuncs = genpath('./dynamicFunctions');
symb = genpath('./symbolicFunctions');
addpath(utilities,symb,dynFuncs);


%% Initialize model parameters
param.m = 1;
k1 = 2;
k2 = 0.2;
param.k = [k1 k2]'; % for setting initial conditions only!
param.c = 0;
n = 3; %state dimension - pos vel inputforce
m = 2; %output dimension - vel force

%% Initialize Filter params
tMin = 0;
tMax = 5;
x_ic = 2.0;     %system position initial condition
v_ic = 0;       %yMeas(1,1);%0; %system velocity initial condition
f_ic  = 0;      % yMeas(1,2);

param.dtKalman = 0.001;      %EKF computation time step (discretisation)

tKalman = tMin:param.dtKalman:tMax;

z0 = [x_ic v_ic f_ic];

param.t = tMin;

model = param;

%% Obtain noisy measurement
dataFolder = './robotData/Measurement/';
load(strcat(dataFolder,'measurement.mat'));
var = variance;

t = tspan;
yMeas =interp1(t,y,tKalman);
z = interp1(t,z,tKalman);
f = interp1(t,f,tKalman);
k1 = interp1(t,k1,tKalman);
k2 = interp1(t,k2,tKalman);
c = interp1(t,c,tKalman);
var_vel = var(1);
var_f = var(2);


%% Function handles

f_func = @forwardDynamics;
h_func = @output;
df_dx_func = @derivativeForwardDynamics;
dh_dx_func = @derivativeOutput;


%% Initialize the covariance matrices,

P = diag([0.001,0.0001,0.001]);    % NxN state covariance matrix
Q = diag([0.0,0.0,0.0025]);   % randomly assumed Process noise covariance - design parameter for the filter performance
% R = diag([var_vel var_f]);
R = diag([50.0 0.25]);%27.5

%% Initialize the state estimate and covariance
xh = z0'; %Initial x0
Ph = P;

Xhat      = zeros(n,length(tKalman))';
Xupdt      = zeros(n,length(tKalman))';

%% Run the filter
for i = 1:length(tKalman)
    if(mod(i,100) == 0)
        fprintf('Count: %d , Time: %f \n', i , toc());
    end
    tic;
    %update step
    [xh,Ph] = ekf_update1(xh,Ph,yMeas(i,:)',dh_dx_func,R,h_func,[],model);
    Xupdt(i,:) = xh;
    Ph = (Ph + Ph')/2;
    P(:,:,i) = Ph;
    
    %predict step
    [xh,Ph] = ekf_predict1(xh,Ph,df_dx_func,Q,f_func,[],model);
    Xhat(i,:) = xh;
    
end




%% Plot the desired behavior, measurement, predicted and corrected estimates
figure(1)
plot(tKalman,z(:,1),'b');
hold on
plot(tKalman,Xupdt(:,1),'g');
title('position vs time');
legend('actual','estimated');

figure(2)
plot(tKalman,yMeas(:,1),'y');
hold on
plot(tKalman,z(:,2),'b');
hold on
plot(tKalman,Xupdt(:,2),'k');
title('velocity vs time');
legend('measurement','actual','estimated');

figure(3)
plot(tKalman,yMeas(:,2),'y');
hold on
plot(tKalman,f','b');
hold on
plot(tKalman,Xupdt(:,3),'k');
title('force vs time');
legend('measurement','actual','estimated');

figure(4)
plot(tKalman,k1,'b');
hold on
plot(tKalman,k2,'r');
hold on
plot(tKalman,c,'g');
title('Stiffness and damping evolution');
legend('k1','k2','c');

