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
param.k1 = 2.0;
param.k2 = 0.2;
param.k = [param.k1 param.k2]';
param.c = 0;
n = 3;
m = 2;
p = 3;


%% Initialize Filter params
tMin = 0;
tMax = 30;
x_ic = 2.0;     %system position initial condition
v_ic = 0;       %yMeas(1,1);%0; %system velocity initial condition
f_ic  = 0;      % yMeas(1,2);
k1_ic = param.k(1);
k2_ic = param.k(2);
c_ic = param.c;


param.dtKalman = 0.001;      %EKF computation time step (discretisation)

tKalman = tMin:param.dtKalman:tMax;

z0 = [x_ic v_ic f_ic];
w0 = [k1_ic k2_ic c_ic];

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

f_func = @forwardDynamicsDualState;
h_func = @outputDualState;
w_func = @paramDynamicsDualState;
df_dx_func = @derivativeForwardDynamicsDualState;
dw_dk_func = @derivativeParamDynamicsDualState;
dh_dx_func = @derivativeOutputDualState;
dc_dw_func = @derivativeParamOutputDualState;


%% Initialize the covariance matrices,

P = diag([0.001,0.0001,0.001]);
P_param = diag([0.02,0.002,0.02]);

Q = diag([0.0,0.0,0.0025]);
Q_param = diag([0.0125,0.0075,0.0007]); %0.25%2.5%3.4975

R = diag([50.0 2.5]);
R_param = R;

%% Initialize the state estimate and covariance
xh = z0'; %Initial x0
wh = w0';
Ph = P;
Pw = P_param;

Xhat      = zeros(n,length(tKalman))';
What      = zeros(p,length(tKalman))';
Xupdt      = zeros(n,length(tKalman))';
Wupdt      = zeros(p,length(tKalman))';

%% Run the filter
for i = 1:length(tKalman)
    if(mod(i,100) == 0)
        fprintf('Count = %d ; Time = %f \n',i,toc());
    end
    %update step
    tic;
     
    
    [xh,Ph] = ekf_update1(xh,Ph,yMeas(i,:)',dh_dx_func,R,h_func,[],model);
    Xupdt(i,:) = xh;
    Ph = (Ph + Ph')/2;
    P(:,:,i) = Ph;
    
    %PARAMETER update : uses a function ekf_updateparam1, modified from the function ekf_update1
    % function [M,P,K,MU,S,LH] = ekf_updateparam1(X,M,P,y,H,C,R,h,V,param)
    % X : state estimate from the state filter
    % M : parameter estimate
    % H : State Measurement jacobian
    % C : Parameter Measurement Jacobian
    [wh, Pw] = ekf_updateparam1(xh,wh,Pw,yMeas(i,:)',dh_dx_func,dc_dw_func,R_param,h_func, [], model);
    Wupdt(i,:) = wh;
    Pw = (Pw + Pw')/2;
    P_param(:,:,i) = Pw;
    model.k(1) = wh(1);
    model.k(2) = wh(2);
    model.c = wh(3);

    
    
    %PARAMETER prediction
    %uses function ekf_predictparam1 accepts state, parameter inputs
    [wh, Pw] = ekf_predictparam1(xh,wh, Pw,dw_dk_func,Q_param,w_func,[], model);
    model.k(1) = wh(1);
    model.k(2) = wh(2);
    model.c = wh(3);
    What(i,:) = wh;
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
plot(tKalman,Wupdt(:,1),'k');
title('k1 vs time');
legend('actual','estimated');

figure(5)
plot(tKalman,k2,'b');
hold on
plot(tKalman,Wupdt(:,2),'k');
title('k2 vs time');
legend('actual','estimated');

figure(6)
plot(tKalman,c,'b');
hold on
plot(tKalman,Wupdt(:,3),'k');
title('c vs time');
legend('actual','estimated');

