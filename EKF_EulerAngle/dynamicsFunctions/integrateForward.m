function dx = integrateForward(t,x,f_B_o,mu_B_o,f_B_c,mu_B_c,model)

v_B = x(1:3);
omega_B = x(4:6);
phi = x(7:9);

G_g = model.G_g;
m = model.m;
I_B = model.I;


k = model.k;
K = diag([k(1),k(2),k(3)]);
c = model.c;
C = diag([c(1),c(2),c(3)]);

phi0 = model.phi0';


dv_B     = -S(omega_B) * v_B + 1/m * f_B_o - 1/m*f_B_c + euler2dcm(phi)*G_g;
domega_B =  I_B \ ((-S(omega_B) * (I_B * omega_B)) + mu_B_o - mu_B_c + (- K'*K*(phi - phi0') - C*omega_B));
dphi     =  Tomega_dphi(phi)\omega_B;

dx = [dv_B;domega_B;dphi];

end

