function dh = derivativeOutputJointState(Xhat, model)


dh = outputsDerivativesJointState(...
    Xhat(1) ,Xhat(2) , Xhat(3) , ...
    model.m, Xhat(4) ,Xhat(5) , Xhat(6));



end
