function dh = derivativeOutput(Xhat, model)


dh = outputsDerivatives(...
    Xhat(1) ,Xhat(2) , Xhat(3) , ...
    model.m, model.k(1),model.k(2),model.c);



end
