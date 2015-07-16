function xn = forwardDynamics( x0, model)
dt = model.dtKalman;

t_u = model.t;


massSpringProcess = @(t,x)processExplicitODE(x0(1),x0(2),model.m,model.k,model.omega,model.Amplitude,t_u);
[~,x] = ode45(massSpringProcess,[0 dt],x0);



xn = x(end,:)';

end