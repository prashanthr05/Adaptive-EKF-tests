function h = outputJointState(Xhat,model)

h = measurementModelJointState(Xhat(1),Xhat(2),Xhat(3),model.m,Xhat(4),Xhat(5),Xhat(6));

end
