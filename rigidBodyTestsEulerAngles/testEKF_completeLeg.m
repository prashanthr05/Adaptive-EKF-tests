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


%% Model Parameters common across experiments
%loads measurement, model, initial conditions, and timescale
dataFolder = './Measurement/';
load(strcat(dataFolder,'measurement.mat'));

model.dtKalman = 0.01; % EKF computation time step (discretisation)

model.measurementPlots = 'noPlots'; % options - 'makePlots' , 'noPlots'
model.filterOutputPlots = 'makePlots';
model.skipSteps = 50; % no of steps to skip for diplaying kalman execution time in loop


%% basic experiment setup
model.filter = 'jekf'; %options 'ekf' 'jekf' 'dekf'


fprintf('Running %s .... \n',model.filter);

if(strcmp(model.filter,'ekf') == 1)
    n       = 21;         % state dimension - (translational vel, rotational vel, w_o, w_c, RPY angle)
    m       = 18;         % output dimension
    
    f_func     = @forwardDynamics;
    df_dx_func = @wrapperForwardDynamicsDerivative;
    dh_dx_func = @wrapperOutputsDerivatives;
    h_func     = @outputFunc;
end

if(strcmp(model.filter,'jekf') == 1)
    n       = 27;         % state dimension - (translational vel, rotational vel, w_o, w_c, RPY angle, params)
    m       = 18;         % output dimension
    
    f_func     = @forwardDynamicsJointState;
    df_dx_func = @wrapperForwardDynamicsDerivativeJointState;
    dh_dx_func = @wrapperOutputsDerivativesJointState;
    h_func     = @outputFuncJointState;
end


if(strcmp(model.filter,'dekf') == 1)
    n = 21;
    p = 6;
    m = 18;
    
    f_func     = @forwardDynamics;
    df_dx_func = @wrapperForwardDynamicsDerivative;
    dh_dx_func = @wrapperOutputsDerivatives;
    h_func     = @outputFunc;
    
    w_func = @paramDynamicsDualState;
    dw_dk_func = @wrapperParamDynamicsDerivative;
    dc_dw_func = @wrapperParamOutputsDerivatives;
end

[kalmanQParams,kalmanRParams,k,c] = setupCovariancesForExperiments(model.filter);
k = model.k';
c = model.c';


%% Kalman Parameters

% setting up process covariances
kalman.a_Q  = kalmanQParams(1);%4.5;
kalman.omega_Q  = kalmanQParams(2);%4.0;%4.75;
kalman.f_Q  = kalmanQParams(3);%0.5;%6.5;
kalman.mu_Q =kalmanQParams(4);%2.5;%6.5;
kalman.phi_Q = kalmanQParams(5);%1.5;%2.50;
kalman.K_Q =kalmanQParams(6);
kalman.C_Q =kalmanQParams(7);

% setting up measurement covariances
kalman.sigma_f = kalmanRParams(1);
kalman.sigma_u = kalmanRParams(2);
kalman.sigma_a = kalmanRParams(3);
kalman.sigma_omega = kalmanRParams(4);
kalman.sigma_skin = kalmanRParams(5);

%% Obtaining measurement data
%    [yMeas,tMeas,model,RData] = realMeasurement_completeLeg(model.dtKalman,model,model.measurementPlots,t_min,t_max,'withSkin',7,'right','right');
%    T = tMeas(end);
%    tKalman = tMeas;
%
t_min = model.t_min;
t_max = model.t_max;

t = linspace(t_min,t_max,(t_max - t_min)/model.dtKalman);
tKalman = t';


for i = 1 : size(y,1)
    yMeas(i,:) = interp1(tDes,y(i,:)',t);
end

for i = 1 : size(sysX,1)
    sys(i,:) = interp1(tDes,sysX(i,:)',t);
end
%% Initializing covariances
if(strcmp(model.filter,'ekf') == 1)
    Q  = diag([kalman.a_Q*ones(3,1);...
        kalman.omega_Q*ones(3,1);...
        kalman.f_Q*ones(3,1);...
        kalman.mu_Q*ones(3,1);...
        kalman.f_Q*ones(3,1);...
        kalman.mu_Q*ones(3,1);...
        kalman.phi_Q*ones(3,1)]);
    
    kalman.P = 15*diag([kalman.a_Q*ones(3,1);kalman.omega_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.phi_Q*ones(3,1)]);
    
    x0 = model.x0;
    
else if(strcmp(model.filter,'jekf') == 1)
        
        Q  = diag([kalman.a_Q*ones(3,1);...
            kalman.omega_Q*ones(3,1);...
            kalman.f_Q*ones(3,1);...
            kalman.mu_Q*ones(3,1);...
            kalman.f_Q*ones(3,1);...
            kalman.mu_Q*ones(3,1);...
            kalman.phi_Q*ones(3,1);...
            kalman.K_Q*ones(3,1);...
            kalman.C_Q*ones(3,1)]);
        
        kalman.P = 15*diag([kalman.a_Q*ones(3,1);kalman.omega_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.phi_Q*ones(3,1);kalman.K_Q*ones(3,1);kalman.C_Q*ones(3,1)]);
        
        x0 = [model.x0;model.w(:,1)];
        %       x0 = [model.x0;0.0;4.5;0.0;0.0;0.5;0.0];
        
    else if(strcmp(model.filter,'dekf') == 1)
            
            Q  = diag([kalman.a_Q*ones(3,1);...
                kalman.omega_Q*ones(3,1);...
                kalman.f_Q*ones(3,1);...
                kalman.mu_Q*ones(3,1);...
                kalman.f_Q*ones(3,1);...
                kalman.mu_Q*ones(3,1);...
                kalman.phi_Q*ones(3,1)]);
            
            Q_param = diag([kalman.K_Q*ones(3,1);...
                kalman.C_Q*ones(3,1)]);
            
            kalman.P  = 15*diag([kalman.a_Q*ones(3,1);...
                                 kalman.omega_Q*ones(3,1);...
                                 kalman.f_Q*ones(3,1);...
                                 kalman.mu_Q*ones(3,1);...
                                 kalman.f_Q*ones(3,1);...
                                 kalman.mu_Q*ones(3,1);...
                                 kalman.phi_Q*ones(3,1)]);
            
            P_param = 15*diag([kalman.K_Q*ones(3,1);...
                               kalman.C_Q*ones(3,1)]);
            R_param = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3)]);
            
            x0 = model.x0;
            w0 = model.w(:,1);
        end
    end
end
forceR = 'true';

if(strcmp(forceR,'true')==1 || ~exist('RData'))
    disp('Assuming an R value');
    
    %     R = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3),kalman.sigma_skin.*ones(1,1)]);
    R = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3)]);
else
    % modify following two lines to check for without skin option
    disp('Using real data covariance matrix');
    RData(18,18) = 35.63;
    R = RData;
end

%% KALMAN FILTER IMPLEMENTATION
%% initialising EKF
% Initializing estimate and update
Ph = kalman.P;
xh        = x0;

if(strcmp(model.filter,'dekf') == 1)
    wh = w0;
    Pw = P_param;
    Xhat      = zeros(n+p,length(tKalman))';
    % Initializing update
    Xupdt = zeros(length(tKalman),n+p);
    P = zeros(n, n,length(tKalman));
    P_p = zeros(p, p,length(tKalman));
else
    Xhat      = zeros(n,length(tKalman))';
    Xupdt = zeros(length(tKalman),n);
    P = zeros(size(Ph,1), size(Ph,2),length(tKalman));
    
end

model.dt = model.dtKalman;


disp('Starting Kalman Filter prediction');
drawnow;

%% EKF execution
for i = 1:length(tKalman)
    tic;
    % Update step
    [xh, Ph] = ekf_update1(xh , Ph, yMeas(:,i), dh_dx_func, R,h_func, [], model);
    
    if(strcmp(model.filter,'dekf') ~= 1)
        Xupdt(i,:) = xh;
    end
    
    Ph = (Ph + Ph')/2;
    xAfterUpdate = xh;
    pAfterUpdate = Ph;
    
    if(strcmp(model.filter,'dekf') == 1)
        %param update
        [wh, Pw] = ekf_updateparam1(xh,wh,Pw,yMeas(:,i),dh_dx_func,dc_dw_func,R_param,h_func, [], model);
        model.w = wh;
        Xupdt(i,:) = [xh;wh];
        Pw = (Pw + Pw')/2';
        
        %param predict
        [wh, Pw] = ekf_predictparam1(xh,wh, Pw,dw_dk_func,Q_param,w_func,[], model);
        model.w = wh;
        
    end
    
    
    % Prediction step
    [xh, Ph] =  ekf_predict1(xh, Ph, df_dx_func, Q, f_func, [], model);
    
    if(strcmp(model.filter,'dekf') ~= 1)
        Xhat(i,:) = xh;
        P(:,:,i)  = Ph;
    else if(strcmp(model.filter,'dekf') == 1)
            
            Xhat(i,:) = [xh;wh];
            P(:,:,i)  = Ph;
            P_p(:,:,i) = Pw;
        end
    end
    
    if(mod(i,model.skipSteps)==0)
        %         fprintf('\n %d Steps processing time : ',model.skipSteps);
        %         disp(toc());
        fprintf('\n Time  : %d \n',i*model.dtKalman);
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

yMeas = yMeas';
sys = sys';
save(strcat(dataBaseFolder,'filteredResult.mat'),'tKalman','yMeas','Xupdt','Xhat','P','sys');
if(strcmp(model.filterOutputPlots,'makePlots') == 1)
    plotAndSaveFigs(dataBaseFolder,plotFigBaseFolder,model.filter);
end





