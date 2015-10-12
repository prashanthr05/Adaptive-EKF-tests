clc
clear 
close all

addpath(genpath('./../utils'));

%% Mass and Inertia parameters of the foot
syms m I_Bxx I_Bxy I_Bxz I_Byy I_Byz I_Bzz real

%% Orientation of the foot expressed in ZYZ euler angle with respect to the inertial frame <O>
syms phi1 phi2 phi3 real 

%% Linear velocity and force components of the body expressed in <B> frame
syms v_Bx v_By v_Bz real
syms f_B_ox f_B_oy f_B_oz real
syms f_B_cx f_B_cy f_B_cz real

%% Angular velocity and moment components of the body expressed in <B> frame
syms omega_Bx omega_By omega_Bz real
syms mu_B_ox mu_B_oy mu_B_oz real
syms mu_B_cx mu_B_cy mu_B_cz real

%% Acceleraton due to gravity expressed in inertial frame <O>
syms G_g1 G_g2 G_g3 real
G_g = [G_g1 G_g2 G_g3]';

%% Compliant Contact Torsional Spring Damper Parameters
syms k_xx k_yy k_zz c_xx c_yy c_zz real
syms phi01 phi02 phi03 real

phi0 = [phi01 phi02 phi03]';

K = diag([k_xx;k_yy;k_zz]);
C = diag([c_xx;c_yy;c_zz]);

w = [k_xx k_yy k_zz c_xx c_yy c_zz]';

%% Parameter Dynamics
dk = zeros(3,1);
dc = zeros(3,1);
dw = [dk;dc];

%% Foot Pose
v_B = [v_Bx v_By v_Bz]';
omega_B = [omega_Bx omega_By omega_Bz]';
phi = [phi1 phi2 phi3]';

%% Wrenches acting on the foot
f_B_o = [f_B_ox f_B_oy f_B_oz]';
mu_B_o = [mu_B_ox mu_B_oy mu_B_oz]';
f_B_c = [f_B_cx f_B_cy f_B_cz]';
mu_B_c = [mu_B_cx mu_B_cy mu_B_cz]';

%% Inertia matrix
I_B = [I_Bxx I_Bxy I_Bxz;...
       I_Bxy I_Byy I_Byz;...
       I_Bxz I_Byz I_Bzz];
   
 dI = [I_Bxx I_Bxy I_Bxz I_Byy I_Byz I_Bzz ]';  
%% Rotation matrix describing orientation of body with respect to the inertial frame <O>
G_R_B = euler2dcm(phi);
   
%% Relation between body angular rate and euler angle rates

   T_B = [cos(phi3)*sin(phi2)  -sin(phi3) 0;...
          sin(phi3)*sin(phi2)   cos(phi3) 0;...
          cos(phi2)                  0        1];
      
%    T_G = [0  -sin(phi1) cos(phi1)*sin(phi2);...
%           0   cos(phi1) sin(phi1)*sin(phi2);...
%           1      0        cos(phi2)];


%% State vector
x_withoutCompliance     = [ v_B;  omega_B;  f_B_o;  mu_B_o; f_B_c;  mu_B_c;  phi];
x_withCompliance     = [ v_B;  omega_B;  f_B_o;  mu_B_o; f_B_c;  mu_B_c;  phi; w];

%% Prediction model
dv_B = -S(omega_B)*v_B + 1/m*f_B_c - 1/m*f_B_o + transpose(G_R_B)*G_g;

domega_B_withCompliance = I_B\(-S(omega_B)*I_B*omega_B + mu_B_c - mu_B_o - K*K'*(phi - phi0) - C*omega_B);
domega_B_withoutCompliance = I_B\(-S(omega_B)*I_B*omega_B + mu_B_c - mu_B_o);

df_B_o = zeros(3,1);
dmu_B_o = zeros(3,1);
df_B_c = zeros(3,1);
dmu_B_c = zeros(3,1);

dphi = T_B\omega_B;
% dphi = T_G\omega_G;

%% Process model Jacobians

f_withCompliance     = [dv_B; domega_B_withCompliance; df_B_o;  dmu_B_o;df_B_c; dmu_B_c; dphi;dw]; %EKF,DEKF
f_withoutCompliance     = [dv_B; domega_B_withoutCompliance; df_B_o;  dmu_B_o;df_B_c; dmu_B_c; dphi]; %EKF


df_dx_withoutCompliance = jacobian(f_withoutCompliance, x_withoutCompliance); 
df_dx_withCompliance = jacobian(f_withCompliance, x_withCompliance);


%% Measurement model 
h_imu = [dv_B - transpose(G_R_B)*G_g; omega_B];
h_fto = [f_B_o; mu_B_o];
h_skin = f_B_cz;
% h_skin = [f_B_cz; mu_B_cy; mu_B_cz];

h_ftc_foot_withoutCompliance = [f_B_c; mu_B_c]; 
h_ftc_foot_withCompliance = [f_B_c; mu_B_c - K*K'*(phi - phi0) - C*omega_B]; 

 

%% Measurement model Jacobians

%Foot No compliance
h_FNC_withSkin = [h_imu ; h_fto; h_ftc_foot_withoutCompliance; h_skin];
h_FNC_withoutSkin = [h_imu ; h_fto; h_ftc_foot_withoutCompliance];

dh_dx_FNC_withSkin = jacobian(h_FNC_withSkin,x_withoutCompliance); 
dh_dx_FNC_withoutSkin = jacobian(h_FNC_withoutSkin,x_withoutCompliance); 

%Foot with compliance
h_FC_withSkin = [h_imu ; h_fto; h_ftc_foot_withCompliance; h_skin];
h_FC_withoutSkin = [h_imu ; h_fto; h_ftc_foot_withCompliance];

dh_dx_FC_withSkin = jacobian(h_FC_withSkin,x_withCompliance); 
dh_dx_FC_withoutSkin = jacobian(h_FC_withoutSkin,x_withCompliance);

model.I  = I_B;

%% Observability test - Foot Without Compliance without skin
% zeroth order lie derivative
dLf0h = dh_dx_FNC_withoutSkin;
% first order lie derivative
Lf1h = dLf0h*f_withoutCompliance;
dLf1h = jacobian(Lf1h,x_withoutCompliance);
% second order lie derivative
Lf2h = dLf1h*f_withoutCompliance;
dLf2h = jacobian(Lf2h,x_withoutCompliance);

O = [dLf0h; dLf1h; dLf2h];
observable_states = rank(O);
unobservable_states = length(x_withoutCompliance) - rank(O);

fprintf('\nObservability of leg without compliance without skin measurement model \n');
fprintf('# Observable States : %d \n', observable_states);
fprintf('# Unobservable States : %d \n', unobservable_states);

%% Observability test - Foot With Compliance
% zeroth order lie derivative
dLf0h = dh_dx_FC_withoutSkin;
% first order lie derivative
Lf1h = dLf0h*f_withCompliance;
dLf1h = jacobian(Lf1h,x_withCompliance);
% second order lie derivative
Lf2h = dLf1h*f_withCompliance;
dLf2h = jacobian(Lf2h,x_withCompliance);

O = [dLf0h; dLf1h; dLf2h];
observable_states = rank(O);
unobservable_states = length(x_withCompliance) - rank(O);

fprintf('\nObservability of leg with compliance without skin measurement model \n');
fprintf('# Observable States : %d \n', observable_states);
fprintf('# Unobservable States : %d \n', unobservable_states);
 
