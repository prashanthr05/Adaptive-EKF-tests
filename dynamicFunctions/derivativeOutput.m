function dh = derivativeOutput(Xhat, model)


dh = outputsDerivatives(...
    Xhat(1) ,Xhat(2) , ...
    model.m, model.k,model.omega,model.Amplitude,model.dtKalman);



end