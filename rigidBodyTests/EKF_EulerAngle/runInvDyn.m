clear
close all
clc

%% adding of paths
utilities    = genpath('./utils');
symb         = genpath('./symbolicFunctions');
ellipses      = genpath('./ellipses');
dynFuncs     = genpath('./dynamicsFunctions');
plotFuncs   = genpath('./plotFunctions');
skinFuncs   = genpath('./skinFunctions');
addpath(utilities, symb, ellipses,dynFuncs,plotFuncs,skinFuncs)

%% basic experiment setup
n       = 21;         % state dimension - (translational vel, rotational vel, w_o, w_c, RPY angle)
m       = 19;         % output dimension

%% Model Parameters common across experiments
setup.dtInvDyn = 0.00001;

setup.t_min = 3.0; % time until which to calibrate
setup.t_max = 7.5; % Max time in dataset until which to filter
setup.measurementPlots = 'noPlots'; % options - 'makePlots' , 'noPlots'
setup.skipSteps = 50; % no of steps to skip for diplaying kalman execution time in loop

[kalmanQParams,kalmanRParams,k,c] = setupCovariancesForExperiments();
setup.k = k; % initial stiffness (defined but not used for all experiments)
setup.c = c; % initial damping (defined but not used for all experiments)
setup.w = [setup.k; setup.c];

%% Obtaining Mass and Inertias
setup.m = 4.9580;

%% obtaining mass-matrix from mexWBIModel (stored in utils)
load('./utils/MIni','M_at_q0');

%% TODO We need to extract foot's instead.
Mjts = M_at_q0(7:end,7:end);            % Leg mass matrix 
%MI_at_thigh = Mjts(14:16,14:16);

%%% checkup if this is for left leg or right leg. Not entirely sure if it will be different
%Tl = [0 1 0; 1 0 0 ; 0 0 -1];
%% Transformation of the mass matrix to COM reference frame
%Ic = Tl*MI_at_thigh*Tl' - S([0 0 -0.18102]')*S([0 0 -0.18102]')';
MI_at_thigh_r = Mjts(20:22,20:22);
 Tl_r = [0 -1 0; 1 0 0 ; 0 0 1];
 Ic = Tl_r*MI_at_thigh_r*Tl_r' - S([0 0 -0.18102]')*S([0 0 -0.18102]')';

% disp(Ic);
setup.I = Ic;
%model.x0 = [zeros(3,1);omega(1,:)';fo(:,1);muo(:,1);fc(:,1);muc(:,1);zeros(3,1)];
setup.phi0 = [0,0.5*pi,0]';
%model.acclSign = -1;

%% world gravity
setup.B0_g = [0;0;9.8];
B_R_G = euler2dcm(setup.phi0);
setup.G_g = -B_R_G'*setup.B0_g;

model = setup;

t_min = model.t_min;
t_max = model.t_max;

[f_B1_t, mu_B1_t, f_B2_t, mu_B2_t, x0] = InverseDynamics(t_max, model, 'true');