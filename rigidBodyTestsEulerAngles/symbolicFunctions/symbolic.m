clear 
close all
clc

%% Utils folder contains the function euler2dcm()
utils =  genpath('./../utils');
addpath(utils);

%% Initialise the symbolic variables
syms I_B f_B_o mu_B_o f_B_c mu_B_c m real
syms I_Bxx I_Byy I_Bzz I_Bxy I_Byz I_Bxz G_g1 G_g2 G_g3  real
syms v_Bx v_By v_Bz real
syms f_B_ox f_B_oy f_B_oz real
syms f_B_cx f_B_cy f_B_cz real
syms omega_Bx omega_By omega_Bz real
syms mu_B_ox  mu_B_oy  mu_B_oz real
syms mu_B_cx  mu_B_cy  mu_B_cz real
syms phi1 phi2 phi3 real
syms K C w real
syms k_xx k_yy k_zz c_xx c_yy c_zz real
syms phi01 phi02 phi03 real % rest position

K = diag([k_xx,k_yy,k_zz]);
C = diag([c_xx,c_yy,c_zz]);
w = [k_xx k_yy k_zz c_xx c_yy c_zz]';

I_B = [...
     I_Bxx I_Bxy I_Bxz;...
     I_Bxy I_Byy I_Byz;...
     I_Bxz I_Byz I_Bzz];
dI = [I_Bxx I_Bxy I_Bxz I_Byy I_Byz I_Bzz ]';

omega_B = [    omega_Bx   omega_By   omega_Bz]';
v_B     = [    v_Bx       v_By       v_Bz]';
f_B_o   = [    f_B_ox     f_B_oy     f_B_oz]';
f_B_c   = [    f_B_cx     f_B_cy     f_B_cz]';
mu_B_o  = [    mu_B_ox    mu_B_oy    mu_B_oz]';
mu_B_c  = [    mu_B_cx    mu_B_cy    mu_B_cz]';
phi     = [    phi1       phi2       phi3]';
B_R_G       = euler2dcm(phi);
G_g = [G_g1;G_g2;G_g3];

phi0    = [    phi01      phi02      phi03]'; %rest orientation used by the spring torque

dv_B     = -S(omega_B) * v_B + 1/m * f_B_o - 1/m*f_B_c + euler2dcm(phi)*G_g;
domega_B =  I_B \ (-S(omega_B) * (I_B * omega_B) + mu_B_o - mu_B_c  + (- K'*K*(phi - phi0) - C*omega_B));
df_B_o   =  [0 0 0]';
dmu_B_o  =  [0 0 0]';
df_B_c   =  [0 0 0]';
dmu_B_c  =  [0 0 0]';
dphi     =  (Tomega_dphi(phi))\omega_B;


h_imu = [(dv_B - B_R_G*G_g); omega_B];
h_fto = [f_B_o; mu_B_o];
h_ftc = [f_B_c; mu_B_c];
h_skin = f_B_cz;


f     = [dv_B; domega_B; df_B_o;  dmu_B_o;df_B_c; dmu_B_c; dphi];
x     = [ v_B;  omega_B;  f_B_o;  mu_B_o; f_B_c;  mu_B_c;  phi];
dw    = zeros(6,1);

f_jointState = [f;dw];

df_dx = jacobian(f, x);
df_dw = jacobian(f,w);
df_dx_jointState = [df_dx df_dw;...
      zeros(size(df_dw')) eye(size(df_dw,2))];

h_withSkin = [h_imu ; h_fto; h_ftc; h_skin];
h_withoutSkin = [h_imu ; h_fto; h_ftc];

dh_dx = jacobian(h_withSkin,x);
dh_dw = jacobian(h_withoutSkin,w);
dh_dx_withoutSkin = jacobian(h_withoutSkin,x);
dh_dx_jointState = [dh_dx_withoutSkin dh_dw];

dp_dw = jacobian(dw,w);
dg_dw = dh_dx_withoutSkin*jacobian(f,w);

model.I  = I_B;

%% Regular EKF

matlabFunction(h_withoutSkin,'file','./../dynamicsFunctions/measurement_withoutSkin','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(f,'file','./../dynamicsFunctions/processImplicitODE','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(df_dx,'file','./dynamicsDerivatives','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(dh_dx_withoutSkin,'file','./outputsDerivatives_withoutSkin','vars',[x; w; phi0; dI; m; G_g]);
% matlabFunction(dh_dx_withoutSkin,'file','./outputsDerivatives_withoutSkin','vars',[x; w; phi0; dI; m; G_g]);

%% JEKF

matlabFunction(f_jointState,'file','./../dynamicsFunctions/processImplicitODEJointState','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(df_dx_jointState,'file','./dynamicsDerivativesJointState','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(dh_dx_jointState,'file','./outputsDerivatives_withoutSkinJointState','vars',[x; w; phi0; dI; m; G_g]);

%% DEKF
matlabFunction(dw,'file','./../dynamicsFunctions/paramExplicitODEDualState','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(dp_dw,'file','./paramdynamicsDerivativesDualState','vars',[x; w; phi0; dI; m; G_g]);
matlabFunction(dg_dw,'file','./paramoutputsDerivativesDualState','vars',[x; w; phi0; dI; m; G_g]);
