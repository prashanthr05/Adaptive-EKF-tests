% Department of Robotics, Brain and Cognitive Sciences
% Istituto Italiano di Tecnologia, 11 September 2014
%
% Original code by: Francesco Nori,
% Modified to include local parametrization by: Jorhabib Eljaik G. and
% Naveen Kuppuswamy
%
% Modified to include constraint moment by: Prashanth Ramadoss
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
dynFuncs     = genpath('./dynamicsFunctions');
plotFuncs   = genpath('./plotFunctions');
skinFuncs   = genpath('./skinFunctions');
addpath(utilities, symb, dynFuncs, plotFuncs, skinFuncs)


%% Model Parameters common across experiments
%setup.dtInvDyn = 0.00001;
setup.dtForDyn = 0.0001; % EKF forward dynamics computation time step
setup.dtKalman = 0.01; % EKF computation time step (discretisation)

setup.urdf = 'icubGazeboSim'; %options 'icubGazeboSim' 'icubGenova03'
% setup.urdf = 'icubGenova03';
% setup.urdf = 'iCubHeidelberg01';


setup.t_min = 3.0; % time until which to calibrate
setup.t_max = 8.0; % Max time in dataset until which to filter
if(strcmp(setup.urdf,'icubGazeboSim') == 1)
    setup.t_min = 3.0; % time until which to calibrate
    setup.t_max = 300.0;%300; % Max time in dataset until which to filter
end

setup.measurementPlots = 'makePlots'; % options - 'makePlots' , 'noPlots'
setup.filterOutputPlots = 'makePlots';
setup.skipSteps = 50; % no of steps to skip for diplaying kalman execution time in loop

expt = 'foot'; %options -'leg', 'foot'
filter = 'ekf';
legChoice = 'left';
skinChoice = 'right';
%% toggle with dataset
if(strcmp(setup.urdf,'icubGazeboSim') == 1)
        dataSet = 'FWD'; % options 'BWD' 'FWD'
        if(strcmp(dataSet,'FWD') == 1)
            expNo = 1; %don't change
        else (strcmp(dataSet,'BWD') == 1)
            expNo = 1; % options 1 2
        end
else
    dataSet = 'old'; % options 'old' 'new' NOTE : for leg.. use only old dataset
    if(strcmp(dataSet,'old') == 1)
        expNo = 7; %don't change
    else (strcmp(dataSet,'new') == 1)
        expNo = 4; % options 1 - firm surface; 2 - onelayermat; 3 - twolayermat; 4 - threelayermat
    end
end


%% Filter EKF
fprintf('\nRunning %s \n-------------\n\n',filter);

clearvars -except filter filterID i setup dataSet expNo expt legChoice skinChoice
close all


if(strcmp(setup.urdf,'icubGazeboSim') == 1)
    [kalmanQParams,kalmanRParams,kIni,cIni] = gazeboCovariancesForExperiments(filter);
else
    [kalmanQParams,kalmanRParams,kIni,cIni] = robotCovariancesForExperiments(filter);
end
setup.kalmanQParams = kalmanQParams;
setup.kalmanRParams = kalmanRParams;
setup.k = kIni; % initial stiffness (defined but not used for all experiments)
setup.c = cIni; % initial damping
setup.w = [kIni;cIni];

setup.filter = 'ekf';

if(strcmp(setup.urdf,'icubGazeboSim') == 1)
    expID = 2;%:2;
    % 1 - w/o Skin w/o compliance
    % 2 - w/0 skin w compliance
    
    measurementSuffix = {'withoutSkin','withoutSkin'};
    processSuffix = {'withoutCompliance','withCompliance'};
    n       = [21,27];         % state dimension - (translational vel, rotational vel, w_o, w_c, RPY angle)
    m       = [18,18];         % output dimension
    p       = [0,0];
else
    expID = 1;
    % 1 - w/o Skin w/o compliance
    % 2 - w/o skin w compliance
    % 3 - w skin w/o compliance
    % 4 - w skin w compliance
    measurementSuffix = {'withoutSkin','withoutSkin','withSkin','withSkin'};
    processSuffix = {'withoutCompliance','withCompliance','withoutCompliance','withCompliance'};
    n       = [21,27,21,27];         % state dimension - (translational vel, rotational vel, w_o, w_c, RPY angle)
    m       = [18,18,19,19 ];         % output dimension
    p       = [0,0,0,0];
end




for j = expID
    
    fprintf('\nProcessing measurement %s with process %s\n-------------\n\n',measurementSuffix{j},processSuffix{j});
    clearvars -except filter filterID i j expID measurementSuffix processSuffix n m p setup dataSet expNo expt legChoice skinChoice
    close all
    
    %% Function handls for process mode, measurement model and their derivatives
    
    f_func     = str2func(strcat('forwardDynamics_',processSuffix{j}));
    df_dx_func = str2func(strcat('wrapperForwardDynamicsDerivative_',processSuffix{j}));
    dh_dx_func = str2func(strcat('wrapperOutputsDerivatives_',expt,'_',processSuffix{j},'_',measurementSuffix{j}));
    h_func = str2func(strcat('output_',expt,'_',processSuffix{j},'_',measurementSuffix{j}));
    
    
    %% Kalman Parameters
    kalmanQParams = setup.kalmanQParams;
    % setting up process covariances
    kalman.a_Q  = kalmanQParams(1);%4.5;
    kalman.omega_Q  = kalmanQParams(2);%4.0;%4.75;
    kalman.f_Q  = kalmanQParams(3);%0.5;%6.5;
    kalman.mu_Q =kalmanQParams(4);%2.5;%6.5;
    kalman.phi_Q = kalmanQParams(5);%1.5;%2.50;
    kalman.k_Q =kalmanQParams(6);
    kalman.c_Q =kalmanQParams(7);
    
    
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
    
    
    if(strcmp(setup.urdf,'icubGazeboSim') == 1)
        [yMeas,tMeas,model] = gazeboMeasurement_comprehensive(model.dtKalman,model,model.measurementPlots,t_min,t_max,measurementSuffix{j},processSuffix{j},expNo,legChoice,skinChoice,dataSet,expt);
    else
        [yMeas,tMeas,model] = robotMeasurement_comprehensive(model.dtKalman,model,model.measurementPlots,t_min,t_max,measurementSuffix{j},processSuffix{j},expNo,legChoice,skinChoice,dataSet,expt);
    end
    T = tMeas(end);
    tKalman = tMeas;
    
    %% Setting up Filter Covariances
    
    %Measurement Noise Covariance
    forceR = 'true';
    
    if(strcmp(forceR,'true')==1 || ~exist('RData'))
        disp('Assuming an R value');
        if(strcmp(measurementSuffix{j},'withoutSkin') == 1)
            R = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3)]);
        else if(strcmp(measurementSuffix{j},'withSkin') == 1)
                R = diag([kalman.sigma_a.*ones(1,3), kalman.sigma_omega.*ones(1,3), kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3),kalman.sigma_f.*ones(1,3), kalman.sigma_u.*ones(1,3),kalman.sigma_skin.*ones(1,1)]);
            end
        end
    else
        
        % modify following two lines to check for without skin option
        disp('Using real data covariance matrix');
        RData(19,19) = 35.63;
        R = RData;
    end
    
    %Process Noise and State Transition Covariance
    switch(processSuffix{j})
        case 'withoutCompliance'
            Q  = diag([kalman.a_Q*ones(3,1);
                kalman.omega_Q*ones(3,1);
                kalman.f_Q*ones(3,1);
                kalman.mu_Q*ones(3,1);
                kalman.f_Q*ones(3,1);
                kalman.mu_Q*ones(3,1);
                kalman.phi_Q*ones(3,1)]);
            kalman.P = 15*diag([kalman.a_Q*ones(3,1);kalman.omega_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.phi_Q*ones(3,1)]);
            x0 = model.x0;
            
        case 'withCompliance'
            Q  = diag([kalman.a_Q*ones(3,1);
                kalman.omega_Q*ones(3,1);
                kalman.f_Q*ones(3,1);
                kalman.mu_Q*ones(3,1);
                kalman.f_Q*ones(3,1);
                kalman.mu_Q*ones(3,1);
                kalman.phi_Q*ones(3,1);
                kalman.k_Q*ones(3,1);
                kalman.c_Q*ones(3,1)]);
            
            kalman.P = 15*diag([kalman.a_Q*ones(3,1);kalman.omega_Q*ones(3,1);kalman.f_Q*ones(3,1);kalman.mu_Q*ones(3,1);kalman.f_Q*ones(3,1);...
                kalman.mu_Q*ones(3,1);kalman.phi_Q*ones(3,1);kalman.k_Q*ones(3,1);kalman.c_Q*ones(3,1)]);
            x0 = model.x0;
    end
    
    
    
    %% KALMAN FILTER IMPLEMENTATION
    %% initialising EKF
    % Initializing estimate and update
    model.dt = model.dtKalman;
    Ph = kalman.P;
    xh        = x0;
    
    
    Xhat      = zeros(n(j),length(tKalman))';
    Xupdt = zeros(length(tKalman),n(j));
    P = zeros(size(Ph,1), size(Ph,2),length(tKalman));
    
    
    
    disp('Starting Kalman Filter prediction');
    drawnow;
    
    %% EKF execution
    for k = 1:length(tKalman)
        tic;
        
        % Update step
        [xh, Ph] = ekf_update1(xh , Ph, yMeas(k,:)', dh_dx_func, R,h_func, [], model);
        
        Xupdt(k,:) = xh;
        Ph = (Ph + Ph')/2;
        xAfterUpdate = xh;
        pAfterUpdate = Ph;
        
        
        % Prediction step
        [xh, Ph] =  ekf_predict1(xh, Ph, df_dx_func, Q, f_func, [], model);
        
        Xhat(k,:) = xh;
        P(:,:,k)  = Ph;
        
        
        if(mod(k,model.skipSteps)==0)
            %                         fprintf('\n %d Steps processing time : ',model.skipSteps);
            %                         disp(toc());
            fprintf('\n Time  : %d \n',t_min + (k*model.dtKalman));
        end
        
    end
    
    
        
    plotFigBaseFolder = sprintf('./plots/%s%s/',processSuffix{j},measurementSuffix{j});
    dataBaseFolder = sprintf('./data/%s%s/',processSuffix{j},measurementSuffix{j});
    
    %% RMSE
    
%     orientationError = model.phi_groundtruth' - Xupdt(:,19:21).*(180/pi);
%     orientationRMSE = rms(orientationError);
% 
%    save(strcat(dataBaseFolder,'errorResultPeriod05.mat'),'tKalman','orientationError','orientationRMSE','Xupdt','model');

%% Plots
    
    if(~exist(dataBaseFolder))
        mkdir(dataBaseFolder);
    end
    
    finalVersion =0;
    if(finalVersion == 1)
        close all;
    end
    
    save(strcat(dataBaseFolder,'filteredResult.mat'),'tKalman','yMeas','Xupdt','Xhat','P','model');
    if(strcmp(model.filterOutputPlots,'makePlots') == 1)
        plotAndSaveFigs(dataBaseFolder,plotFigBaseFolder,processSuffix{j});
    end
    
    if(j < 4)
%         fprintf('\nPress any key to continue to next experiment\n');
%         pause;
    end
    
    
end



