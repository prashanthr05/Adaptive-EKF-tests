function [model,tMax,leg_ft,foot_ft,inertial,transforms] = gazeboSensorParams_comprehensive(whichLeg,whichSkin,numberOfExperiment,t_max,dtKalman,model,dataSet,expt)
%% function returns the model parameters for treating the complete leg as a single rigid body.
%The function for now only returns the values corresponding to right leg

%% Setting datasets

% expPath     = strcat('./robotData/gazebo',dataSet,'Tipping',num2str(numberOfExperiment),'/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/gazeboFWDTipping1/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/leftFootRaise/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/gazeboStanding1/dumper/icubGazeboSim/');
expPath     = strcat('./robotData/gazeboBWDTipping1/dumper/icubGazeboSim/');

% expPath     = strcat('./robotData/fallTest/dumper/icubGazeboSim/');

% expPath     = strcat('./robotData/LeftFootPitchRotatePeriod0.01/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/LeftFootPitchRotatePeriod0.1/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/LeftFootPitchRotatePeriod1/dumper/icubGazeboSim/');

% expPath     = strcat('./robotData/lfp05/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/lfp01/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/lfp005/dumper/icubGazeboSim/');



% expPath     = strcat('./robotData/LeftFootRollRotatePeriod0.01/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/LeftFootRollRotatePeriod0.1/dumper/icubGazeboSim/');
% expPath     = strcat('./robotData/LeftFootRollRotatePeriod1/dumper/icubGazeboSim/');

% expPath     = strcat('./robotData/KneeJointRotate/dumper/icubGazeboSim/');

leg_choice  = whichLeg;


% % Leg F/T analog sensor
leg_ft_data   = importdata(strcat(expPath,leg_choice,'_leg/analog:o/data.log'));
% Foot F/T analog sensor
foot_ft_data  = importdata(strcat(expPath,leg_choice,'_foot/analog:o/data.log'));
% Inertial sensor attached to the foot
inertial_data = importdata(strcat(expPath,leg_choice,'_foot_IMU/data.log'));


leg_ft.t = leg_ft_data(:,2)-leg_ft_data(1,2);
leg_ft.idx = leg_ft_data(:,1) - leg_ft_data(1,1);
leg_ft.data(:,3:8) = leg_ft_data(:,3:8);
leg_ft.f = leg_ft.data(:,3:5);
leg_ft.mu = leg_ft.data(:,6:8);

%%
foot_ft.t = foot_ft_data(:,2)-foot_ft_data(1,2);
foot_ft.idx = foot_ft_data(:,1) - foot_ft_data(1,1);
foot_ft.data(:,3:8) = foot_ft_data(:,3:8);
foot_ft.f = foot_ft.data(:,3:5);
foot_ft.mu = foot_ft.data(:,6:8);


inertial.t = inertial_data(:,2)-inertial_data(1,2);
inertial.idx = inertial_data(:,1) - inertial_data(1,1);
inertial.data = inertial_data(:,3:end);


ttmp = linspace(foot_ft.t(1),foot_ft.t(end),length(inertial.t));
inertial.data = interp1(ttmp,inertial.data,foot_ft.t);
inertial.t = foot_ft.t;
inertial.idx = foot_ft.idx;

%% Transformation between IMU sensor frame and body coordinate frame located at CoM
com_R_imu = eye(3);
transforms.B_R_imu = com_R_imu;

%% TODO - include leg transformation to foot
if(strcmp(expt,'foot')==1)
%     [model,transforms,leg_ft] = dynComputations_foot(model,inertial,leg_ft,transforms,leg_choice,expPath,t_max,dtKalman);
    [model,transforms,leg_ft] = dynComputations_foot_nolegFT(model,inertial,leg_ft,transforms,leg_choice,expPath,t_max,dtKalman);
    
    ttmp1 = linspace(foot_ft.t(1),foot_ft.t(end),length(leg_ft.t));
    leg_ft.f = interp1(ttmp1,leg_ft.f,foot_ft.t);
    leg_ft.mu = interp1(ttmp1,leg_ft.mu,foot_ft.t);
    leg_ft.t = foot_ft.t;
    leg_ft.idx = foot_ft.idx;
    
else if(strcmp(expt,'leg')==1)
        [model,transforms] = dynComputations_leg(model,inertial,transforms,leg_choice,expPath);
    end
end


%% world gravity
model.B0_g = [0;0;9.8];
G_R_B = euler2dcm(model.phi0);
model.G_g = G_R_B*model.B0_g;


%%
tMax = min([foot_ft.t(end),leg_ft.t(end),inertial.t(end),t_max]);
%tMax = min([foot_ft.t(end),inertial.t(end),t_max]);
% tMax = min([leg_ft.t(end),foot_ft.t(end),skin.t(end),inertial.t(end),t_max]);


end

