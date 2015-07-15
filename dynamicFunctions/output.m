function h = output(Xhat,model)

h = measurementModel(Xhat(1),Xhat(2),model.m,model.k,model.omega,model.Amplitude,model.dtKalman);

end
