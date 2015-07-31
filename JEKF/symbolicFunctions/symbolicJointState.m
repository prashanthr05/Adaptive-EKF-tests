clear 
close all
clc

syms x v X real
syms m k k1 k2 c real
syms F real

X = [x v F]'; %position velocity inputforce
k = [k1 k2]';
w = [k; c];

dx = X(2);
dv = (-k(1)*X(1) - k(2)*X(1)*X(1)*X(1) - c*dx + F)/m;
dF = 0;

dX = [dx dv dF]';

dK =  zeros(2,1); %randomwalk
dc = zeros(1,1);  %randomwalk

dW = [dK; dc];

f = dX;
h = [v F]';

df_dX = jacobian(f,X);
df_dw = jacobian(f,w);
dh_dX = jacobian(h,X);
dh_dw = jacobian(h,w);

f_jointState = [dX; dW];
h_jointState = h;

df_jointState = [df_dX df_dw;...
  zeros(size(df_dw')) eye(size(df_dw,2))];
dh_jointState = [dh_dX dh_dw];

matlabFunction(h_jointState,'file','./../dynamicFunctions/measurementModelJointState','vars',[X;m;k;c]);
matlabFunction(f_jointState,'file','./../dynamicFunctions/processExplicitODEJointState','vars',[X;m;k;c]);
matlabFunction(df_jointState,'file','./dynamicsDerivativesJointState','vars',[X;m;k;c]);
matlabFunction(dh_jointState,'file','./outputsDerivativesJointState','vars',[X;m;k;c]);
