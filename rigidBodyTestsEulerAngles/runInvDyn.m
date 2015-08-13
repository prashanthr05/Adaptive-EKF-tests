% This script is used to generate the forces using an inverse dynamics function
% The forces are later used in a forward dynamics function to verify the
% evolution of equations of motion. Helps in verifying the correctness of
% our dynamic model.

% All the model parameters are initialised in this script. This calls the
% invDyn(..) function which is used to set the trajectories of orientation
% and velocity of rigid body. rigidBodyInvDyn(..) computes the net torques
% acting on the body for the given acce3lerations and velocities.


clear
close all
clc

%% adding of paths
utilities    = genpath('./utils');
dynFunctions = genpath('./dynamicsFunctions');
addpath(utilities,dynFunctions);

%% Model Parameters common across experiments
setup.dtInvDyn = 1e-4;
setup.t_topple = 6.5;
setup.t_min = 0.0; % time until which to calibrate
setup.t_max = 12.0; % Max time in dataset until which to filter

setup.invDynPlots = 'makePlots'; %makePlots noPlots
setup.measurementPlots = 'noPlots';
%% Obtaining Mass and Inertias
setup.m = 4.9580;

%% obtaining mass-matrix from mexWBIModel (stored in utils)
load('./utils/MIni','M_at_q0');
Mjts = M_at_q0(7:end,7:end);            % Leg mass matrix
MI_at_thigh_r = Mjts(20:22,20:22);
Tl_r = [0 -1 0; 1 0 0 ; 0 0 1];
Ic = Tl_r*MI_at_thigh_r*Tl_r' - S([0 0 -0.18102]')*S([0 0 -0.18102]')';

setup.I = Ic;

%% world gravity
setup.phi0 = [0;0;0]; % just for definition ---> set properly in invDyn
setup.B0_g = [0;0;9.8]; % just for definition ---> set properly in invDyn
setup.B_R_G = euler2dcm(setup.phi0); % just for definition ---> set properly in invDyn
setup.G_g = setup.B_R_G'*setup.B0_g; % just for definition ---> set properly in invDyn

%% Setting up stiffness and damping
setup.k = [0.0,5,0.0]; %just for definition ---> set properly in invDyn
setup.c = [0.0,5.0,0.0]; %just for definition ---> set properly in invDyn
setup.w = [setup.k';setup.c'];

model = setup;

t_min = model.t_min;
t_max = model.t_max;

tDes = 0:0.001:t_max;

%% Orientation and Velocity trajectory response - should be differentiable functions
rotTraj = 'simReality'; %options - 'step' 'ramp' 'sine' 'simReality'; rotTraj is set only to the Pitch angle  
velTraj = 'simReality'; %options - 'step' 'ramp' 'd_ellipse' 'simReality'; velTraj is set to only along negative x-axis and  positive z-axis
param = 'varying'; %options - 'varying' 'not'
%% Obtaining Noisy measurementes
[f_B1_t, mu_B1_t, f_B2_t, mu_B2_t, omega_B, ddPos, dPos, Phi, B_g, w, x0,model] = invDyn(t_max, model,tDes,rotTraj,velTraj,param);

model.x0 = x0;
model.w = w(tDes);
a_B = ddPos(tDes);
v_B = dPos(tDes);
omega_B = omega_B(tDes);
f_B1 = f_B1_t(tDes);
mu_B1 = mu_B1_t(tDes);
f_B2 = f_B2_t(tDes);
mu_B2 = mu_B2_t(tDes);
Phi_B = Phi(tDes);
B_g = B_g(tDes);

% Output from the inverse dynamics for comparison with forward dynamics and filter results
sysX = [v_B; omega_B; f_B1; mu_B1; f_B2; mu_B2; Phi_B; model.w];  %repmat(model.w,1,size(v_B,2))

% Force measurement with known variance white noise
variance_f = 0.001;
wn_f = sqrt(variance_f).*randn(size(f_B1)); % white noise with a given variance
f1 = f_B1 + wn_f;
f2 = f_B2 + wn_f;

% Torque measurement with known variance white noise
variance_mu = 0.001;
wn_mu = sqrt(variance_mu).*randn(size(mu_B1)); % white noise with a given variance
mu_1 = mu_B1 + wn_mu;
mu_2 = mu_B2 + wn_mu;

% Angular Velocity measurement with known variance white noise
variance_omega = 0.0001;
wn_omega = sqrt(variance_omega).*randn(size(omega_B)); % white noise with a given variance
omega = omega_B + wn_omega;

% Normal Force due to the Skin
fcZ = f2(3,:);

% Measurement Vector
% y = [a_B - B_g; omega; f1; mu_1; f2; mu_2; fcZ]; %acc omega fo muo fc muc fcz
 y = [a_B - B_g; omega; f1; mu_1; f2; mu_2];
%% Saving the data : used by filter and forward Dynamics
dataFolder = sprintf('./Measurement/');
save(strcat(dataFolder,'measurement.mat'),'y','model','tDes','sysX');


%% Plot InvDyn results
if(strcmp(model.invDynPlots,'makePlots') == 1)
    figure(1);
    subplot(2,1,1);
    plot(tDes,f_B1); hold on;
    legend('fx','fy','fz');
    title('Inverse Dynamics - Forces due to body')
    subplot(2,1,2);
    plot(tDes,f_B2); hold on;
    legend('fx','fy','fz');
    title('Inverse Dynamics - Forces due to contact')
    
    figure(2);
    subplot(2,1,1);
    plot(tDes,mu_B1); hold on;
    legend('mux','muy','muz');
    title('Inverse Dynamics - Moments due to body')
    subplot(2,1,2);
    plot(tDes,mu_B2); hold on;
    legend('mux','muy','muz');
    title('Inverse Dynamics - Moments due to contact')
    
    
    figure(3);
    subplot(2,1,1);
    plot(tDes,a_B); hold on;
    legend('ax','ay','az');
    title('Acc');
    subplot(2,1,2);
    plot(tDes,omega_B); hold on;
    legend('\omega_{x}','\omega_{y}','\omega_{z}');
    title('Angular vel');
    
    figure(4)
    subplot(2,1,1)
    plot(tDes,Phi_B); hold on;
    legend('phi1','phi2','phi3');
    title('Orientation');
    subplot(2,1,2)
    plot(tDes,v_B); hold on;
    legend('vx','vy','vz');
    title('Lin vel');
end

%% Plot simulated measurements
if(strcmp(model.measurementPlots,'makePlots') == 1)
    figure(5);
    subplot(2,1,1);
    plot(tDes,f1); hold on;
    legend('fx','fy','fz');
    title('Measurement - Forces due to body')
    subplot(2,1,2);
    plot(tDes,f2); hold on;
    legend('fx','fy','fz');
    title('Measurement - Forces due to contact')
    
    figure(6);
    subplot(2,1,1);
    plot(tDes,mu_1); hold on;
    legend('mux','muy','muz');
    title('Measuremnt - Moments due to body')
    subplot(2,1,2);
    plot(tDes,mu_2); hold on;
    legend('mux','muy','muz');
    title('Measurement - Moments due to contact')
    
    
    figure(7);
    subplot(2,1,1);
    plot(tDes,a_B - B_g); hold on;
    legend('ax','ay','az');
    title('Acc measured');
    subplot(2,1,2);
    plot(tDes,omega); hold on;
    legend('\omega_{x}','\omega_{y}','\omega_{z}');
    title('Angular velocity - measured');
    
    figure(8);
    plot(tDes,fcZ); hold on;
    title('fcz Skin');
    
end
