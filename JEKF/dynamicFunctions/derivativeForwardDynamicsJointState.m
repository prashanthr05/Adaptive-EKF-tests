function A = derivativeForwardDynamicsJointState(Xhat, model)

dt = model.dtKalman;

df = dynamicsDerivativesJointState(...
    Xhat(1) ,Xhat(2),Xhat(3) , ...
    model.m, Xhat(4),Xhat(5),Xhat(6));

A = df.*dt + eye(size(df));

end
