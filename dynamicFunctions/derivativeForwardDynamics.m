function A = derivativeForwardDynamics(Xhat, model)

dt = model.dtKalman;

df = dynamicsDerivatives(...
    Xhat(1) ,Xhat(2) , ...
    model.m, model.k,model.omega,model.Amplitude,dt);

A = df.*dt + eye(size(df));

end