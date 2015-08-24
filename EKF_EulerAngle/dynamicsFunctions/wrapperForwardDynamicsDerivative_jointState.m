function A = wrapperForwardDynamicsDerivative_jointState(Xhat, model)

dt = model.dt;



df = dynamicsDerivatives_jointState(...
    Xhat(1) ,Xhat(2) ,Xhat(3),...
    Xhat(4) ,Xhat(5) ,Xhat(6),...
    Xhat(7) ,Xhat(8) ,Xhat(9),...
    Xhat(10),Xhat(11),Xhat(12),...
    Xhat(13),Xhat(14),Xhat(15),...
    Xhat(16),Xhat(17),Xhat(18),...
    Xhat(19),Xhat(20),Xhat(21),...
    Xhat(22),Xhat(23),Xhat(24),Xhat(25),Xhat(26),Xhat(27),...
    model.phi0(1),model.phi0(2),model.phi0(3),...
    model.I(1,1), model.I(1,2),model.I(1,3),model.I(2,2),model.I(2,3),model.I(3,3),model.m,model.G_g(1),model.G_g(2),model.G_g(3));





A = df.*dt + eye(size(df));
