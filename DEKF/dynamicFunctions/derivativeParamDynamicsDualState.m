function A = derivativeParamDynamicsDualState(Xhat, model)

dt = model.dtKalman;

 df = paramdynamicsDerivativesDualState(...
     Xhat(1) ,Xhat(2), Xhat(3) , ...
     model.m, model.k(1),model.k(2),model.c);

A = df.*dt + eye(size(df));

end
