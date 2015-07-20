clear 
close all
clc


syms x v X real
syms m k k1 k2 real
syms omega Amplitude F real
syms t real

X = [x v]';
k = [k1 k2]';
dx = X(2);
dv = (-k(1)*X(1) - sign(X(1))*k(2)*X(1)*X(1) + Amplitude*sin(omega*t))/m;
% dv = (-k(1)*X(1) - k(2)*X(1)*X(1)*X(1) + Amplitude*sin(omega*t))/m;
dX = [dx dv]';

f = dX;
h = v;
df_dX = jacobian(f,X);
dh_dX = jacobian(h,X);


matlabFunction(df_dX,'file','./dynamicsDerivativesNonLinear','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(dh_dX,'file','./outputsDerivativesNonLinear','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(f,'file','./../dynamicFunctions/processExplicitODENonLinear','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(h,'file','./../dynamicFunctions/measurementModelNonLinear','vars',[X;m;k;omega;Amplitude;t]);