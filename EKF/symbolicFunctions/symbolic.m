clear 
close all
clc

syms x v X real
syms m k k1 k2 c real
syms F real

% X = [x v]';
X = [x v F]'; %position velocity inputforce
k = [k1 k2]';
dx = X(2);
% dv = (-k(1)*X(1) - sign(X(1))*k(2)*X(1)*X(1) + F)/m;
dF = 0;
dv = (-k(1)*X(1) - k(2)*X(1)*X(1)*X(1) - c*dx + F)/m;
dX = [dx dv dF]';
% dX = [dx dv]';
dK =  zeros(2,1);

f = dX;


h = [v F]';
df_dX = jacobian(f,X);
df_dw = jacobian(f,k);
dh_dX = jacobian(h,X);
dh_dw = jacobian(h,k);

matlabFunction(df_dX,'file','./dynamicsDerivatives','vars',[X;m;k;c]);
matlabFunction(dh_dX,'file','./outputsDerivatives','vars',[X;m;k;c]);
matlabFunction(f,'file','./../dynamicFunctions/processExplicitODE','vars',[X;m;k;c]);
matlabFunction(h,'file','./../dynamicFunctions/measurementModel','vars',[X;m;k;c]);
