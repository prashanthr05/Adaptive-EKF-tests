function dg_dw = paramoutputsDerivativesDualState(x,v,F,m,k1,k2,c)
%PARAMOUTPUTSDERIVATIVESDUALSTATE
%    DG_DW = PARAMOUTPUTSDERIVATIVESDUALSTATE(X,V,F,M,K1,K2,C)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    31-Jul-2015 04:00:37

t2 = 1.0./m;
t3 = x.^2;
dg_dw = reshape([-t2.*x,0.0,-t2.*t3.*x,0.0,-t2.*v,0.0],[2, 3]);
