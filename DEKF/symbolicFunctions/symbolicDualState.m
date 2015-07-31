clear 
close all
clc


syms x v X real
syms m k k1 k2 c real
syms  F real

 X = [x v F]'; %position velocity inputforce

k = [k1 k2]';
w = [k; c];

dx = X(2);
dv = (-k(1)*X(1) - k(2)*X(1)*X(1)*X(1) - c*dx + F)/m;
dF = 0;

dX = [dx dv dF]';


dw =  zeros(3,1);
f = dX;
h = [v F]';


df_dX = jacobian(f,X);
dp_dw = jacobian(dw,w);
dh_dX = jacobian(h,X);
dg_dw = dh_dX*jacobian(f,w);


matlabFunction(df_dX,'file','./dynamicsDerivativesDualState','vars',[X;m;k;c]);
matlabFunction(dh_dX,'file','./outputsDerivativesDualState','vars',[X;m;k;c]);
matlabFunction(f,'file','./../dynamicFunctions/processExplicitODEDualState','vars',[X;m;k;c]);
matlabFunction(h,'file','./../dynamicFunctions/measurementModelDualState','vars',[X;m;k;c]);
matlabFunction(dw,'file','./../dynamicFunctions/paramExplicitODEDualState','vars',[X;m;k;c]);
matlabFunction(dp_dw,'file','./paramdynamicsDerivativesDualState','vars',[X;m;k;c]); 
matlabFunction(dg_dw,'file','./paramoutputsDerivativesDualState','vars',[X;m;k;c]); 
