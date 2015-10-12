function [model,transforms] = dynComputations_leg(model,inertial,transforms,leg_choice,expPath);

%Defining object of DynamicsComputations
dynComp = iDynTree.DynamicsComputations();

%% Load robot urdf
%Loading urdf of corresponding robot
loadModel = dynComp.loadRobotModelFromFile(strcat(model.urdf,'.urdf'),'urdf'); % Loading iCubGenova03 model.urdf
regCheck = dynComp.isValid(); %checking if regressor generator is correctly configured
dofInternal = dynComp.getNrOfDegreesOfFreedom();

%% Initialize link and frame indices
%Getting indices for base link frame - *_foot, leg FT frame - *_hip_3,
%world frame - locatd at *_sole at time t = 0, when robot is in upright
%standing position. (* = 'l' or 'r')
if(strcmp(leg_choice,'left')== 1)
    footCheck = dynComp.setFloatingBase('l_foot'); %Setting floating base to the l_foot
    base_link = dynComp.getFloatingBase();
    idx = dynComp.getFrameIndex(base_link);
    Nb = dynComp.getFrameIndex('l_hip_3');
%     world = dynComp.getFrameIndex('l_sole');
    world = dynComp.getFrameIndex('l_foot_dh_frame');
else if(strcmp(leg_choice,'right')== 1)
        footCheck = dynComp.setFloatingBase('r_foot'); %Setting floating base to the l_foot
        base_link = dynComp.getFloatingBase();
        idx = dynComp.getFrameIndex(base_link);
        Nb = dynComp.getFrameIndex('r_hip_3');
%         world = dynComp.getFrameIndex('r_sole');
        world = dynComp.getFrameIndex('r_foot_dh_frame');
    end
end


%% Set robot state at time t = 0
%Obtain joint states from torso, right leg and left leg
torso_state = importdata(strcat(expPath,'torso/state:o/data.log'));
torso.timeStamp = torso_state(1,2) - torso_state(1,2) ;
torso.q = torso_state(1,3:end).*(pi/180);

right_leg_state = importdata(strcat(expPath,'right_leg/state:o/data.log'));
right_leg.timeStamp = right_leg_state(1,2) - right_leg_state(1,2) ;
right_leg.q = right_leg_state(1,3:end).*(pi/180);

left_leg_state = importdata(strcat(expPath,'left_leg/state:o/data.log'));
left_leg.timeStamp = left_leg_state(1,2) - left_leg_state(1,2) ;
left_leg.q = left_leg_state(1,3:end).*(pi/180);

%Obtaining acceleration and angular velocity from IMU
com_R_imu = transforms.B_R_imu;
a = (com_R_imu*inertial.data(1,4:6)')';
omega = (com_R_imu*inertial.data(1,7:9)')';

q = iDynTree.VectorDynSize();
q_dot = iDynTree.VectorDynSize();
q_dotdot = iDynTree.VectorDynSize();
base_vel = iDynTree.Twist();
base_acc = iDynTree.ClassicalAcc();
world_T_base = iDynTree.Transform();
world_gravity = iDynTree.SpatialAcc();

jointpos = zeros(dofInternal,1);
jointpos(1:6,1) = left_leg.q(1,1:6);
jointpos(7:12,1) = right_leg.q(1,1:6);
jointpos(15,1) = torso.q(1,1);
jointpos(14,1) = torso.q(1,2);
jointpos(13,1) = torso.q(1,3);
q.fromMatlab(jointpos);

jointvel = zeros(dofInternal,1);
q_dot.fromMatlab(jointvel);
jointacc = zeros(dofInternal,1);
q_dotdot.fromMatlab(jointacc);

lin_vel = zeros(3,1);
ang_vel = omega';
base_vel.fromMatlab([lin_vel;ang_vel]);

lin_acc = a';
ang_acc = zeros(3,1);
base_acc.fromMatlab([lin_acc;ang_acc]);

state_set = dynComp.setRobotState(q,q_dot,q_dotdot,world_T_base,base_vel,base_acc,world_gravity);

%% Leg mass
% Getting the whole leg mass and inertia matrix considering it as a single rigid body
%neglecting foot mass(foot masss = 0.5935).

% Ifoot = dynComp.getLinkInertia(base_link);
% mfoot = Ifoot.getMass(); % mass of the foot
% Icfoot = Ifoot.getRotationalInertiaWrtCenterOfMass();

% foot_p_B = Ifoot.getCenterOfMass();
% B_p_foot = -foot_p_B.toMatlab();

Ileg = iDynTree.SpatialInertia();
for j = idx - 1 : -1 : Nb  
    
    Il = dynComp.getLinkInertia(j);
    foot_X_Nb = dynComp.getRelativeTransform(idx,j);
  
    Ij = foot_X_Nb*Il;
    Ileg = Ileg + Ij;
end


mleg = Ileg.getMass();
I = Ileg.getRotationalInertiaWrtCenterOfMass();
Ic = I.toMatlab();

%% Transforms
%Obtain the center of mass of the leg body wrt to legFT and footFt and
%corresponding adjoint Transforms
foot_p_comleg = Ileg.getCenterOfMass();
ankle_p_B = foot_p_comleg.toMatlab();



legFT_T_foot = dynComp.getRelativeTransform(Nb,idx);
p_comleg = iDynTree.Position(ankle_p_B(1),ankle_p_B(2),ankle_p_B(3));
upperleg_p_comleg = legFT_T_foot*p_comleg;
leg_p_B = upperleg_p_comleg.toMatlab();


model.m = mleg;
Ic = eye(3)*Ic*transpose(eye(3)) - S(ankle_p_B)*S(ankle_p_B)';
model.I = Ic;


transforms.B_adjT_ankle = [eye(3) zeros(3) ; -eye(3)*S(ankle_p_B) eye(3) ];
transforms.B_adjT_leg = [eye(3) zeros(3) ; -eye(3)*S(leg_p_B) eye(3) ];


%% world2body rotation
%Getting world to body rotation, considering foot FT frame 
% and world frame to be located at l_sole but symmetrical in rotation wrt body frame
%To obtain initial body configuration wrt world frame at time t = 0. 
if(state_set)
    
    %getting transform between l_sole and footFT
    footFT_X_sole = dynComp.getRelativeTransform(idx,world);
    footFT_R_sole = footFT_X_sole.getRotation();
   
    % Adding a custom frame at the origin of l_sole but orientation
    % symmetrical to legFT frame and call it world frame/global frame
    % This is done so that we have a fixed inertial frame in which the
    % initial conditions can be set and the gravity is known.
    footFT_p_sole = footFT_X_sole.getPosition(); 
%     footFT_R_sole = iDynTree.Rotation(0,0,1,0,1,0,-1,0,0);
    
    footFT_R_world = iDynTree.Transform(footFT_R_sole,footFT_p_sole);
    
    
    B_R_G = footFT_R_world.getRotation().toMatlab();
    G_R_B = B_R_G';
%     model.phi0 = dcm2euler(G_R_B);
    model.phi0 = [0; pi/2; 0];
end



end
