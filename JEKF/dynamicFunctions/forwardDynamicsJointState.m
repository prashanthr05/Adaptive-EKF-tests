function xn = forwardDynamicsJointState( x0, model)
dt = model.dtKalman;
odeSettings = odeset('InitialStep', 1e-8, 'MaxStep',1e-4);
massSpringProcess = @(t,x)processExplicitODEJointState(x(1),x(2),x(3),model.m,x(4),x(5),x(6));
[~,x] = ode45(massSpringProcess,[0 dt],x0,odeSettings);



xn = x(end,:)';

end
