% function foot_ft = footComputations(leg_ft,t,inertial,transforms,numberOfExpriment)
clc %to be removed
clear %tbr
close all; %tbr
t_min = 0; %tbr
dtKalman = 0.01; %tbr
numberOfExperment = 7; %to be removed

dynComp = iDynTree.DynamicsComputations();
%% Adding robot model iCubGenova03
util = genpath('./utils'); %to be removd 
addpath(util); %to be removd 

modelFolder = './utils/'; 

loadModel = dynComp.loadRobotModelFromFile(strcat(modelFolder,'model.urdf'),'urdf');
regCheck = dynComp.isValid(); %checking if regressor generator is correctly configured
l_footCheck = dynComp.setFloatingBase('l_foot');
base_link = dynComp.getFloatingBase();
dofInternal = dynComp.getNrOfDegreesOfFreedom;
Ic = dynComp.getLinkInertia(base_link);

% initialSamplesToDrop = 100;
samplesToTrim = 50;

expPath     = ['./robotData/backwardTipping/dumperTippingSetup07/icub/']; %to be removed 
%% Obtaining joint velocities and accelerations
torso_state = importdata(strcat(expPath,'torso/state:o/data.log'));
torso.timeStamp = torso_state(:,2) - torso_state(1,2) ;
torso.q = torso_state(:,3:end).*(pi/180);
torso.dq = zeros(size(torso.timeStamp));
torso.ddq = zeros(size(torso.timeStamp));

right_leg_state = importdata(strcat(expPath,'right_leg/state:o/data.log'));
right_leg.timeStamp = right_leg_state(:,2) - right_leg_state(1,2) ;
right_leg.q = right_leg_state(:,3:end).*(pi/180);
right_leg.dq = zeros(size(right_leg.timeStamp));
right_leg.ddq = zeros(size(right_leg.timeStamp));


left_leg_state = importdata(strcat(expPath,'left_leg/state:o/data.log'));
left_leg.timeStamp = left_leg_state(:,2) - left_leg_state(1,2) ;
left_leg.q = left_leg_state(:,3:end).*(pi/180);
left_leg.dq = zeros(size(left_leg.timeStamp));
left_leg.ddq = zeros(size(left_leg.timeStamp));


torso = smoothAndEstimateVelAcc(torso);
right_leg = smoothAndEstimateVelAcc(right_leg);
left_leg = smoothAndEstimateVelAcc(left_leg);

% torsoState = trimDataset(torsoSmooth,samplesToTrim)'

%% Obtaining angular acceleration
inertial_data = importdata(strcat(expPath,'inertial/data.log')); %to be removed
inertial.t = inertial_data(:,2)-inertial_data(1,2); %to be removed
inertial.data = inertial_data(:,3:end); %to be removed

% TODO - express in Foot COM reference frame by properly computing IMU
% rotation
load('IMUOffset_New.mat','com_R_imu'); %to be removed


omega.q = inertial.data(:,7:9);
omega.q = (com_R_imu*omega.q')';
omega.timeStamp = inertial.t;
omega.dq = zeros(size(omega.timeStamp));
omega.ddq = zeros(size(omega.timeStamp));
omega = smoothAndEstimateVelAcc(omega);


acc.q = inertial.data(:,4:6);
acc.q = (com_R_imu*acc.q')';
acc.timeStamp = inertial.t;
acc.dq = zeros(size(acc.timeStamp));
acc.ddq = zeros(size(acc.timeStamp));
acc = smoothAndEstimateVelAcc(acc);


torso = trimDataset(torso,samplesToTrim);
left_leg = trimDataset(left_leg,samplesToTrim);
right_leg = trimDataset(right_leg,samplesToTrim);
omega = trimDataset(omega,samplesToTrim);
acc = trimDataset(acc,samplesToTrim);

t = min([length(acc.timeStamp),length(omega.timeStamp),length(torso.timeStamp),length(right_leg.timeStamp),length(left_leg.timeStamp)]);


q = iDynTree.VectorDynSize();
q_dot = iDynTree.VectorDynSize();
q_dotdot = iDynTree.VectorDynSize();
base_vel = iDynTree.Twist();
base_acc = iDynTree.ClassicalAcc();
world_T_base = iDynTree.Transform();
world_gravity = iDynTree.SpatialAcc(); %by default set to 0, because base acceleration includes gravity components
lin_vel = zeros(3,1);

for i = 1 %: t

    jointpos = zeros(dofInternal,1);
    jointpos(1:6,1) = left_leg.q(i,1:6);
    jointpos(7:12,1) = right_leg.q(i,1:6);
    jointpos(15,1) = torso.q(i,1);
    jointpos(14,1) = torso.q(i,2);
    jointpos(13,1) = torso.q(i,3);
    q.fromMatlab(jointpos);    
       
    jointvel = zeros(dofInternal,1);
    jointvel(1:6,1) = left_leg.dq(i,1:6);
    jointvel(7:12,1) = right_leg.dq(i,1:6);
    jointvel(15,1) = torso.dq(i,1);
    jointvel(14,1) = torso.dq(i,2);
    jointvel(13,1) = torso.dq(i,3);
    q_dot.fromMatlab(jointvel);
        
    jointacc = zeros(dofInternal,1);
    jointacc(1:6,1) = left_leg.ddq(i,1:6);
    jointacc(7:12,1) = right_leg.ddq(i,1:6);
    jointacc(15,1) = torso.ddq(i,1);
    jointacc(14,1) = torso.ddq(i,2);
    jointacc(13,1) = torso.ddq(i,3);
    q_dotdot.fromMatlab(jointacc);
    
    ang_vel = omega.q(i,1:3)';
    base_vel.fromMatlab([lin_vel;ang_vel]);
        
    lin_acc = acc.q(i,1:3)';
    ang_acc = omega.dq(i,1:3)';
    base_acc.fromMatlab([lin_acc;ang_acc]);
    
      
    state_check(i) = dynComp.setRobotState(q,q_dot,q_dotdot,world_T_base,base_vel,base_acc,world_gravity);
    if(state_check(i) == 1)
      idx = dynComp.getFrameIndex(base_link);
      Nb = dynComp.getFrameIndex('l_hip_2');
      
      footFT_X_legFT = dynComp.getRelativeTransform(idx,Nb);
      X_FT = footFT_X_legFT.asAdjointTransformWrench();
%       f_foot(i,1:6) =  footFT_X_legFT * f_leg;
      
      for j = idx - 1 : -1 : Nb
          footFT_X_Nb = dynComp.getRelativeTransform(idx,j);
          X_Nb = footFT_X_Nb.asAdjointTransformWrench();
          
          Il = dynComp.getLinkInertia(j);
          m = Il.getMass();
          I = Il.getRotationalInertiaWrtCenterOfMass();
          I = I.toMatlab();
          
          a = dynComp.getFrameProperSpatialAcceleration(j);
          a_link = a.toMatlab();
          v_dot = a_link(1:3);
          omega_dot = a_link(4:6);
          
          v = dynComp.getFrameTwist(j);
          v_link = v.toMatlab();
          v = v_link(1:3);
          omega = v_link(4:6);
         
          h_dot = [m*eye(3)*v_dot + S(omega)*m*v; I*omega_dot + S(omega)*I*omega]; 
          hdot = iDynTree.VectorDynSize();
          hdot.fromMatlab(h_dot);
%          f_foot(i,1:6) = footFT_X_Nb * hdot;
%           hdot = I*a + S(v.toMatlab())*I*v;
           f_foot(i,1:6) = f_foot(i,1:6) - footFT_X_Nb * hdot;
      end
      
    end
        
    
end

