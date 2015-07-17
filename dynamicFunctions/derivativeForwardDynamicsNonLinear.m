function A = derivativeForwardDynamicsNonLinear(Xhat, model)

dt = model.dtKalman;

df = dynamicsDerivativesNonLinear(...
    Xhat(1) ,Xhat(2) , ...
    model.m, model.k(1),model.k(2),model.omega,model.Amplitude,dt);

A = df.*dt + eye(size(df));

end
