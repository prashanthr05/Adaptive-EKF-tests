clear all
close all
clc

syms x v X real
syms m k real
syms omega Amplitude F real
syms t real

X = [x v]';

dx = X(2);
dv = (-k*X(1) + Amplitude*sin(omega*t))/m;

dX = [dx dv]';

f = dX;
h = v;
df_dX = jacobian(f,X);
dh_dX = jacobian(h,X);


matlabFunction(df_dX,'file','./dynamicsDerivatives','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(dh_dX,'file','./outputsDerivatives','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(f,'file','./../dynamicFunctions/processExplicitODE','vars',[X;m;k;omega;Amplitude;t]);
matlabFunction(h,'file','./../dynamicFunctions/measurementModel','vars',[X;m;k;omega;Amplitude;t]);