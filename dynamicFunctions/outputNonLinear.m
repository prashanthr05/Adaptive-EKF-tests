function h = outputNonLinear(Xhat,model)

h = measurementModelNonLinear(Xhat(1),Xhat(2),model.m,model.k(1),model.k(2),model.omega,model.Amplitude,model.dtKalman);

end