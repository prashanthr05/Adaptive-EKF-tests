function [f_Bo_t, mu_Bo_t, f_Bc_t, mu_Bc_t, omega_B, ddPos,B_g,w] = rigidBodyInvDyn(Phi,dPhi, dPos,ddPos,tc,tDes,model, coeff_f, coeff_mu)
% RIGIDBODYINVDYN This function computes inverse dynamics for a rigid 
% body system with two external forces, assuming null linear velocity and 
% a desired orientation trajectory parametrized in ZYZ  Euler angles (roll,
% pitch and yaw).
%     
% It takes as inputs: 
% Phi   : Orientation as a time series vector
% Pos  : Position of the CoM in inertial frame expressed in body frame as
% time series vector
% model : Model structure containing mass (m), gravity (g), Inertia
%         matrix (I) and the integration time interval (dt).
% c_f12 : Forces ratio
% c_mu12: Torques ratio
%
% Authors: Jorhabib Eljaik,Naveen Kuppuswamy
% Edited by - Prashanth Ramadoss
% Istituto Italiano di Tecnologia
% Department of Robotics, Brain and Cognitive Sciences (RBCS)


    %% Retrieving model parameters
    m   = model.m;
    G_g   = model.G_g;
    I_B = model.I;
    k = model.k;
%     K = diag([k(1),k(2),k(3)]);
    c = model.c;
%     C = diag([c(1),c(2),c(3)]);
    
    Phi0 = model.phi0';
    dt  = model.dtInvDyn;
  
    %% Angular velocity and its derivative given a time-changing orientation 
    % expressed in ZYZ Euler angles
    %Phi      = Phi';

    % Derivative of the time-varying orientation
        T_t      = Tomega_dphi(Phi,'t');
    
    omega_B = zeros(size(Phi));
    domega_B = zeros(size(Phi));
    
        
    % Angular velocity for the desired time-varying orientation
    for i=1:length(dPhi)
         omega_B(i,:) = (T_t{i}*dPhi(i,:)')';
    end

     
    % Derivative of the angular velocity
    domega_B(1:end-1,1:3) = diff(omega_B)/dt;
    domega_B(end,:) = domega_B(end-1,:);

 
    %% Other quantities
    R = euler2dcm(Phi,'t');

    %% External wrenches to be applied
    % Forces
    %f_B1_t = [];
    %f_B2_t = [];
    f_t =zeros(3,length(tc));
    for i=1:length(tc)
       f_t(:,i)   = (m*ddPos(i,:)' + S(omega_B(i,:)')*m*dPos(i,:)' - m*R{i}*G_g);
       B_g(:,i) = R{i}*G_g;
    end
     
    f_Bo_tt = coeff_f*f_t';
    f_Bc_tt = (coeff_f - 1)*f_t';
    
    % Torques
    %mu_B1_t = [];
    %mu_B2_t = [];
    mu_t = zeros(3,length(tc));
    for i=1:length(tc)
        
        mu_t(:,i)  = I_B*domega_B(i,:)' + S(omega_B(i,:)')*I_B*omega_B(i,:)' - (- diag(k(i,:))'*diag(k(i,:))*(Phi(i,:)' - Phi0') - diag(c(i,:))*omega_B(i,:)'); 
          
    end
    mu_Bo_tt = coeff_mu*mu_t';
    mu_Bc_tt = (coeff_mu - 1)*mu_t';
      
   
    f_Bo_t = @(tDes)interp1(tc,f_Bo_tt,tDes)';
    f_Bc_t = @(tDes)interp1(tc,f_Bc_tt,tDes)';
    mu_Bo_t = @(tDes)interp1(tc,mu_Bo_tt,tDes)';
    mu_Bc_t = @(tDes)interp1(tc,mu_Bc_tt,tDes)';
    ddPos = @(tDes)interp1(tc,ddPos,tDes)';
    omega_B = @(tDes)interp1(tc,omega_B,tDes)';
    B_g = @(tDes)interp1(tc,B_g',tDes)';
    w = @(tDes)interp1(tc,model.w,tDes)';
end