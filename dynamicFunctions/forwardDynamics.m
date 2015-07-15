function xn = forwardDynamics( x0, model)
dt = model.dtKalman;

massSpringProcess = @(t,x)processExplicitODE(x0(1),x0(2),model.m,model.k,model.omega,model.Amplitude,t);

[~,x] = ode45(massSpringProcess,[0 dt],x0);
xn = x(end,:)';

end