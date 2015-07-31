function f = processExplicitODE(x,v,F,m,k1,k2,c)
%PROCESSEXPLICITODE
%    F = PROCESSEXPLICITODE(X,V,F,M,K1,K2,C)

%    This function was generated by the Symbolic Math Toolbox version 6.2.
%    30-Jul-2015 17:04:09

t2 = x.^2;
f = [v;-(-F+c.*v+k1.*x+k2.*t2.*x)./m;0.0];
