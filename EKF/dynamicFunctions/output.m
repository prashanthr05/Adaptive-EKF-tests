function h = output(Xhat,model)

h = measurementModel(Xhat(1),Xhat(2),Xhat(3),model.m,model.k(1),model.k(2),model.c);

end
