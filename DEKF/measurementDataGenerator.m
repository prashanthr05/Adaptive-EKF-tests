%measurementDataGenerator
close all;
clear;
clc

%% Initialize model parameters
param.m = 1;
param.k1 = 2;
param.k2 = 0.2;
param.k = [param.k1 param.k2]';
param.c = 1;

param.k1_omega = 0.157;
param.k2_omega = 0.157;
param.k1_C = param.k1;      %Constant bias
param.k2_C = param.k2;      %Constant bias
param.k1_P = 0.6;          %varying part
param.k2_P = 0.06;         %varying part


%% Input Force
param.omega = 1.57;%3.14; % if omega is a non-zero value, it's a forced system
param.Amplitude = 0.5;%2; %Amplitude for forced input

%% Basic params
param.dtForDyn = 0.0001; %EKF computation time step (discretisation)

tMin = 0;
tMax = 30;
x_ic = 2; %system position initial condition
v_ic = 0; %system velocity initial condition
tspan = tMin:param.dtForDyn:tMax;
z0 = [x_ic v_ic];

%% Simulating systm for understanding the true behavior
springmass = @(t,z) [z(2);(-k_varying(param.k1_C, param.k1_P, param.k1_omega,t,tMax)*z(1) - k_varying(param.k2_C, param.k2_P, param.k2_omega,t,tMax)*z(1)*z(1)*z(1) - param.c*z(2) + force(t,param.Amplitude,param.omega,tMax))/param.m];
[t,z] = ode45(springmass,tspan,z0);

%% Evolution of stiffess and damping
for i = 1:length(tspan)
    k1(i) = k_varying(param.k1_C, param.k1_P, param.k1_omega,tspan(i),tMax);
    k2(i) = k_varying(param.k2_C, param.k2_P, param.k2_omega,tspan(i),tMax);
end

c = param.c*ones(size(tspan));

%% Getting velocity measurment
v_offset = 5;
variance_vel = 0.005;
wn_vel = sqrt(variance_vel).*randn(size(z(:,2))); % white noise with a given variance
vMeas = z(:,2) + wn_vel+v_offset;

%% Getting Force measurement
for i = 1:length(tspan)
    f(i) = force(tspan(i),param.Amplitude,param.omega,tMax);
end

variance_f = 0.001;
if(param.Amplitude == 0)
    wn_f = 0;
else
    wn_f = sqrt(variance_f).*randn(size(z(:,2))); % white noise with a given variance
end
fMeas = f' + wn_f;

%% Saving the data in 'measurement.mat'
y = [vMeas fMeas];
variance = [variance_vel variance_f];
t = tspan;
dataFolder = sprintf('./robotData/Measurement/');
save(strcat(dataFolder,'measurement.mat'),'z','f','k1','k2','c','tspan','y','variance');
