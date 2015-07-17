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


%% Initialize model parameters = at the moment run symbolicNonlinear.m to set system to NonLinear and symbolic.m to linear
param.system = 'nonlinear';
param.m = 1;
if(strcmp(param.system,'linear'))
    param.k = 1;
else if(strcmp(param.system,'nonlinear'))
        k1 = 5;
        k2 = 6;
        param.k = [k1 k2]';
    end
end
param.omega = 0; % if omega is a non-zero value, it's a forced system
param.Amplitude = 0; %Amplitude for forced input

param.dtKalman = 0.01; %EKF computation time step (discretisation)


n = 2;
m = 1;

%% Function handles
if(strcmp(param.system,'linear'))
    f_func = @forwardDynamics;
    h_func = @output;
    df_dx_func = @derivativeForwardDynamics;
    dh_dx_func = @derivativeOutput;
else if(strcmp(param.system,'nonlinear'))    
        f_func = @forwardDynamicsNonLinear;
        h_func = @outputNonLinear;
        df_dx_func = @derivativeForwardDynamicsNonLinear;
        dh_dx_func = @derivativeOutputNonLinear;
    end
end

%% Obtain ground truth - To comment if no need to simulate the system
tMin = 0;
tMax = 5;
x_ic = 2; %system position initial condition
v_ic = 0; %system velocity initial condition
tspan = tMin:param.dtKalman:tMax;
z0 = [x_ic v_ic];

param.t = tMin;
if(strcmp(param.system,'linear'))
        springmass = @(t,z) [z(2);(-param.k*z(1) + param.Amplitude*sin(param.omega*t))/param.m];    
        %used as a reference for process model for the EKF
        [t,z] = ode45(springmass,tspan,z0);

else if(strcmp(param.system,'nonlinear'))
        springmass = @(t,z) [z(2);(-param.k(1)*z(1) -sign(z(1))*param.k(2)*z(1)*z(1) + param.Amplitude*sin(param.omega*t))/param.m];    
        %used as a reference for process model for the EKF
        [t,z] = ode45(springmass,tspan,z0);

    end
end



variance = 5.220;
wn = sqrt(variance).*randn(size(z(:,2))); % white noise with a given variance
y = z(:,2) + wn;

% dataFolder = sprintf('./robotData/velocityMeasurement/');
% save(strcat(dataFolder,'measurement2.mat'),'t','y','variance');


%% Obtain noisy measurement
% dataFolder = './robotData/velocityMeasurement/';
% load(strcat(dataFolder,'measurement2.mat'));
var = variance;
tKalman = tspan;
yMeas =y;

model = param;
%% Initialize the covariance matrices
P = 0.01*ones(n,n);%zeros(n,n);%2*eye(n); % NxN state covariance matrix
Q = zeros(n,n);%zeros(n,n);%25*eye(n);    % randomly assumed Process noise covariance
R = var;

%% Initialize the state estimate and covariance
xh = z0'; %Initial x0
Ph = P;

%% Run the filter
for i = 1:length(tKalman)
    %update step
    [xh,Ph] = ekf_update1(xh,Ph,yMeas(i),dh_dx_func,R,h_func,[],model);
    Xupdt(i,:) = xh;
    Ph = (Ph + Ph')/2;
    
    model.t = model.t + model.dtKalman;
    %predict step
    [xh,Ph] = ekf_predict1(xh,Ph,df_dx_func,Q,f_func,[],model);
    Xhat(i,:) = xh;
    P(:,:,i) = Ph;
end




%% Plot the desired behavior, measurement, predicted and corrected estimates
figure(1)
plot(t,z(:,1),'b');
hold on
pause;
plot(t,Xhat(:,1),'r');
hold on
pause;
plot(t,Xupdt(:,1),'g');
% legend('predicted','corrected');
title('position vs time');
legend('actual','predicted','corrected');

figure(2)
plot(t,yMeas,'y');
hold on
pause;
plot(t,z(:,2),'b');
hold on
pause;
plot(t,Xhat(:,2),'r');
hold on
pause;
plot(t,Xupdt(:,2),'k');
% legend('measurement','predicted','corrected');
title('velocity vs time');
legend('actual','measurement','predicted','corrected');
