% Department of Robotics, Brain and Cognitive Sciences
% Istituto Italiano di Tecnologia, 11 September 2014
%
% Original code by: Francesco Nori,
% Modified to include local parametrization by: Jorhabib Eljaik G. and
% Naveen Kuppuswamy
%
% This piece of code simulates the problem of estimating dynamic quantities
% for a single rigid body with distributed force/torque measurements and
% distributed gyro/accelerometers measurements. The motion is governed by
% the following differential equation with a local parametrization of the
% orientation in ZYZ Euler angles.
%
% m    dv^B    + S(omega^B) (m       v^B) = f^B_1  + ... + f^B_n + mg^B
%
% I^B domega^B + S(omega^B) (I^B omega^B) = mu^B_1 + ... + mu^B_n
%
%                                  dphi   = T_phi^-1 * omega
%
% where we defined the following quantities:
%
% I^B    : inertia in the body reference frame
% m      : mass of the rigid body
% f^B_i  : i-th force expressed in the body reference frame
% mu^B_i : i-th torque expressed in the body reference frame
% omega^B: angular velocity expressed in the body reference frame
% v^B    : linear velocity in the body reference frame
% phi    : ZYZ Euler angles representing orientation.
% T_phi  : Transformation matrix between omeg%diff(torsoVelDash)./repmat(diff(t(1:end-1),1,1),1,3);


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
setup.dtForDyn = 0.0001; % EKF forward dynamics computation time step
setup.dtKalman = 0.01; % EKF computation time step (discretisation)

setup.t_min = 3.0; % time until which to calibrate
setup.t_max = 7.5; % Max time in dataset until which to filter
setup.measurementPlots = 'noPlots'; % options - 'makePlots' , 'noPlots'
setup.filterOutputPlots = 'makePlots';
setup.skipSteps = 50; % no of steps to skip for diplaying kalman execution time in loop

[kalmanQParams,kalmanRParams,k,c] = setupCovariancesForExperiments();
setup.kalmanQParams = kalmanQParams;
setup.kalmanRParams = kalmanRParams;
setup.k = k; % initial stiffness (defined but not used for all experiments)
setup.c = c; % initial damping (defined but not used for all experiments)
setup.w = [setup.k; setup.c];


%% Measurement model and its derivative
f_func     = @forwardDynamics;
df_dx_func = @wrapperForwardDynamicsDerivative;
dh_dx_func = @wrapperOutputsDerivatives_withSkin;
h_func     = @outputFunc;

%% Kalman Parameters
kalmanQParams = setup.kalmanQParams;

% setting up process covariances
kalman.a_Q  = kalmanQParams(1);%4.5;
kalman.omega_Q  = kalmanQParams(2);%4.0;%4.75;
kalman.f_Q  = kalmanQParams(3);%0.5;%6.5;
kalman.mu_Q =kalmanQParams(4);%2.5;%6.5;
kalman.phi_Q = kalmanQParams(5);%1.5;%2.50;
kalman.K_Q =kalmanQParams(6);

kalmanRParams = setup.kalmanRParams;

% setting up measurement covariances
kalman.sigma_f = kalmanRParams(1);
kalman.sigma_u = kalmanRParams(2);
kalman.sigma_a = kalmanRParams(3);
kalman.sigma_omega = kalmanRParams(4);
kalman.sigma_skin = kalmanRParams(5);

model = setup;
t_min = model.t_min;
t_max = model.t_max;

%% Obtaining measurement data
   [yMeas,tMeas,model,RData] = realMeasurement_completeLeg(model.dtKalman,model,model.measurementPlots,t_min,t_max,'withSkin',7,'right','right');
   T = tMeas(end);
   tKalman = tMeas;
   
%% Inverse dynamics
% [f_B1_t, mu_B1_t, f_B2_t, mu_B2_t, x0] = InverseDynamics(t_max, model, 'true');


%% Initializing covariances

Q  = diag([kalman.a_Q*ones(3,1);
    kalman.omega_Q*ones(3,1);
    kalman.f_Q*ones(6,1);
    kalman.mu_Q*ones(6,1);
    kalman.phi_Q*ones(3,1)]);

kalman.P = 15*diag([kalman.a_Q*ones(3,1);kalman.omega_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.phi_Q*ones(3,1)]);

forceR = 'true';

if(strcmp(forceR,'true')==1 || ~exist('RData'))
    disp('Assuming an R value');
    
    R = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3),kalman.sigma_skin.*ones(1,1)]);
    
else
    % modify following two lines to check for without skin option
    disp('Using real data covariance matrix');
    RData(19,19) = 35.63;
    R = RData;
end

%% KALMAN FILTER IMPLEMENTATION
%% initialising EKF
% Initializing estimate
Ph = kalman.P;
xh        = model.x0;
Xhat      = zeros(n,length(tKalman))';

% Initializing update
model.dt = model.dtKalman;
Xupdt = zeros(length(tKalman),n);
P = zeros(size(Ph,1), size(Ph,2),length(tKalman));


disp('Starting Kalman Filter prediction');
drawnow;

%% EKF execution
for i = 1:length(tKalman)
    tic;
    
    % Update step
    [xh, Ph] = ekf_update1(xh , Ph, yMeas(i,:)', dh_dx_func, R,h_func, [], model);
    Xupdt(i,:) = xh;
    xAfterUpdate = xh;
    pAfterUpdate = Ph;
    
    % Prediction step
    [xh, Ph] =  ekf_predict1(xh, Ph, df_dx_func, Q, f_func, [], model);
    
    Xhat(i,:) = xh;
    P(:,:,i)  = Ph;
    if(mod(i,setup.skipSteps)==0)
        fprintf('%d Steps processing time :',setup.skipSteps);
        disp(toc());
    end
    
end    
    
    plotFigBaseFolder = sprintf('./plots/');
    dataBaseFolder = sprintf('./data/');
    
    if(~exist(dataBaseFolder))
        mkdir(dataBaseFolder);
    end
    
    finalVersion =0;
    if(finalVersion == 1)
        close all;
    end
    
    save(strcat(dataBaseFolder,'filteredResult.mat'),'tKalman','yMeas','Xupdt','Xhat','P');
    if(strcmp(setup.filterOutputPlots,'makePlots') == 1)
        plotAndSaveFigs(dataBaseFolder,plotFigBaseFolder);
    end
    
    



