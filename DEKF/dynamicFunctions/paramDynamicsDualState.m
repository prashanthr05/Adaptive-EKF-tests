function wn = paramDynamicsDualState(x0, model)
dt = model.dtKalman;
odeSettings      = odeset('InitialStep', 1e-8, 'MaxStep', 1e-4);


massSpringParams = @(t,w)paramExplicitODEDualState(x0(1),x0(2),x0(3),model.m,w(1),w(2),w(3));
[~,w] = ode45(massSpringParams,[0 dt],[model.k(1) model.k(2) model.c],odeSettings);


wn = w(end,:)';

end
