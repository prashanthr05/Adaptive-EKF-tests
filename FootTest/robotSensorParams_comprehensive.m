function [model,tMax,leg_ft,foot_ft,inertial,skin,transforms] = robotSensorParams_comprehensive(whichLeg,whichSkin,numberOfExperiment,t_min,t_max,dtKalman,model,dataSet,expt,measurementType)
%% function returns the model parameters for treating the complete leg as a single rigid body.
%The function for now only returns the values corresponding to right leg

%% Setting datasets

expPath     = strcat('./robotData/backwardTipping/dumperTippingSetup07/icub/');

leg_choice  = whichLeg;
skin_choice = whichSkin;

%% Offsets for the BWD tipping old dataset
    %     Left Leg
    offsets.FT2 = [-4.11735 111.148 -131.104 -13.423 1.36924 -1.04304]; %old dataset
    %     Left Foot
    offsets.FT3 = [-52.0296 -5.64247 -18.0923 -0.516052 9.91345 0.717567]; %old dataset
    %     Right Leg
    offsets.FT4 = [-5.65084 -10.5031 -24.6174 1.95779 2.61084 0.0036913]; %old dataset
    %     Right Foot
    offsets.FT5 = [-17.4992 -0.910681 18.0613 -0.824831 3.32657 -0.252851]; %old dataset


left_leg_ft_offset = offsets.FT2;
right_leg_ft_offset = offsets.FT4;
left_foot_ft_offset = offsets.FT3;
right_foot_ft_offset = offsets.FT5;


if(strcmp(leg_choice,'left')==1)
    leg_ft_offset = left_leg_ft_offset;
    foot_ft_offset = left_foot_ft_offset;
else
    leg_ft_offset = right_leg_ft_offset;
    foot_ft_offset = right_foot_ft_offset;
end
 
% % Leg F/T analog sensor
leg_ft_data   = importdata(strcat(expPath,leg_choice,'_leg/analog:o/data.log'));
% Foot F/T analog sensor
foot_ft_data  = importdata(strcat(expPath,leg_choice,'_foot/analog:o/data.log'));
% Inertial sensor attached to the foot
% inertial_data = importdata(strcat(expPath,leg_choice,'_foot_IMU/data.log'));
inertial_data = importdata(strcat(expPath,'inertial/data.log'));

% Skin Sensor under the foot
	skin_data     = importdata(strcat(expPath,'skin/',skin_choice,'_foot/data.log'));
	skin.t = skin_data(:,2) - skin_data(1,2);
	skin.idx = skin_data(:,1) - skin_data(1,1);
	skin.data = skin_data(:,3:end);


leg_ft.t = leg_ft_data(:,2)-leg_ft_data(1,2);
leg_ft.idx = leg_ft_data(:,1) - leg_ft_data(1,1);
leg_ft.data(:,3:8) = leg_ft_data(:,3:8) - repmat(leg_ft_offset,size(leg_ft_data,1),1);;
leg_ft.f = leg_ft.data(:,3:5);
leg_ft.mu = leg_ft.data(:,6:8);

%%
foot_ft.t = foot_ft_data(:,2)-foot_ft_data(1,2);
foot_ft.idx = foot_ft_data(:,1) - foot_ft_data(1,1);
foot_ft.data(:,3:8) = foot_ft_data(:,3:8) - repmat(foot_ft_offset,size(foot_ft_data,1),1);
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
tCalib = linspace(0,t_min,(t_min)/dtKalman);
a_omega_pre_calib = interp1(inertial.t,inertial.data,tCalib);
a_mean = mean(a_omega_pre_calib(:,4:6)',2);

com_R_imu = computeOptimalIMURotation(dataSet,numberOfExperiment,a_mean,'noPlots','on');
transforms.B_R_imu = com_R_imu;

%% include leg transformation to foot
if(strcmp(expt,'foot')==1)
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
if(strcmp(measurementType,'withSkin')==1)
	tMax = min([leg_ft.t(end),foot_ft.t(end),skin.t(end),inertial.t(end),t_max]);
end


end

