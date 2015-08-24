function [f_B1_t, mu_B1_t, f_B2_t, mu_B2_t, omega_B, ddPos, dPos, Phi B_g,w, x0,model] = InverseDynamics(t_max, model,tDes,rotTraj,velTraj,param)
% This code computes  inverse dynamics for a rigid body system assuming some linear
% velocity and a desired orientation trajectory parametrized in ZYZ
% Euler angles (roll, pitch and yaw).
% Author: Naveen Kuppuswamy
% Edited by : Prashanth Ramadoss
% Istituto Italiano di Tecnologia
% Department of Robotics, Brain and Cognitive Sciences (RBCS)

%% Desired orientation trajectory
dtInvDyn = model.dtInvDyn;
tc = 0:dtInvDyn:t_max;

tsteady = 0:dtInvDyn:model.t_topple;
tdynamic = model.t_topple+dtInvDyn:dtInvDyn:t_max;
tb = dtInvDyn : dtInvDyn : tdynamic(end) - model.t_topple;
%% Setting Orientation trajectory
% % Step Trajectory response
if(strcmp(rotTraj,'step')== 1)
    alpha = zeros(length(tc),1);
    upsilon =  repmat(0.5*pi,length(tc),1)';
    psi = zeros(length(tc),1);
    
    dalpha = zeros(length(tc),1);
    dupsilon = zeros(length(tc),1)';
    dpsi = zeros(length(tc),1);
    
    % Ramp Trajectory Response
else if(strcmp(rotTraj,'ramp')== 1)
        a = 0.35*pi;
        alpha = zeros(length(tc),1);
        upsilon =  0.5*pi + (a/t_max)*tc; % coefficient of the second ramp part must be set in a way to avoid singularities
        psi = zeros(length(tc),1);
        
        dalpha = zeros(length(tc),1);
        dupsilon = repmat((a/t_max),length(tc),1)';
        dpsi = zeros(length(tc),1);
        
        % % sinusoidal trajectory
    else if(strcmp(rotTraj,'sine')== 1)
            Phi_theta =(0.5*pi) /t_max; % singularities lie at angles 0 and pi
            alpha = zeros(length(tc),1);
            upsilon =  0.05+0.9*sin(Phi_theta*tc);
            psi = zeros(length(tc),1);
            
            dalpha = zeros(length(tc),1);
            dupsilon = 0.9*Phi_theta*cos(Phi_theta*tc);
            dpsi = zeros(length(tc),1);
            
            
        else if(strcmp(rotTraj,'simReality')== 1)
                alpha = zeros(length(tc),1);
                
                upsilon1 = repmat(0.5*pi,length(tsteady),1)';
                upsilon2 = upsilon1(end) + (pi/16)*tb.^2;%sin((pi/4)*(tb));
                upsilon = [upsilon1 upsilon2];
                
                psi = zeros(length(tc),1);
                
                dalpha = zeros(length(tc),1);
                
                dupsilon1 = zeros(length(tsteady),1)';
                dupsilon2 = 2*(pi/16)*(tb);
                dupsilon = [dupsilon1 dupsilon2];
                
                dpsi = zeros(length(tc),1);
                
            else
                disp('\n ERR: proper trajecory not set \n')
            end
        end
    end
end
Phi     = [    alpha      upsilon'    psi     ];
dPhi = [ dalpha  dupsilon'  dpsi ];


%% Setting velocity trajectory
% %Assuming step trajectory response
if(strcmp(velTraj,'step')== 1)
    dPos = [-repmat(0.5*pi,length(tc),1),zeros(length(tc),1),repmat(0.5*pi,length(tc),1)];
    ddPos = [-zeros(length(tc),1),zeros(length(tc),1),zeros(length(tc),1)];
    
    % %Assuming ramp trajectory response
else if(strcmp(velTraj,'ramp')== 1)
        dPos = [-(0.5*pi*tc)',zeros(length(tc),1),(0.5*pi*tc)'];
        ddPos = [-repmat(0.5*pi,length(tc),1),zeros(length(tc),1),repmat(0.5*pi,length(tc),1)];
        
        
        % % %Assuming foot position traces an ellipse of (semi major axis a, semi minor
        % % % %axis b).
    else if(strcmp(velTraj,'d_ellipse')== 1)
            ell_a = 3.5;
            ell_b = 1.0;
            ell_th = ((pi/t_max) - ((0.5*pi/t_max)*tc)  )';
            d_th = -(0.5*pi)/t_max;
            
            
            %Pos = [ell_a*cos(ell_th) ; zeros(length(tc),1); ell_b*sin(ell_th)];
            dPos =  [-ell_a*d_th*sin(ell_th)     , zeros(length(tc),1) ,  ell_b*d_th*cos(ell_th)];
            ddPos = [-ell_a*(d_th^2)*cos(ell_th), zeros(length(tc),1) , -ell_b*(d_th^2)*sin(ell_th)];
        else if(strcmp(velTraj,'simReality')== 1)
                dPos1 = [-repmat(0.0,length(tsteady),1),zeros(length(tsteady),1),repmat(0.0,length(tsteady),1)];
                dPos2 = [ -0.25*tb.^2',zeros(length(tb),1),0.5*tb.^2'];
                dPos = [dPos1;dPos2];
                
                ddPos1 = [zeros(length(tsteady),1),zeros(length(tsteady),1),zeros(length(tsteady),1)];
                ddPos2 = [-2*0.25*tb',zeros(length(tb),1),2*0.5*tb'];
                ddPos = [ddPos1; ddPos2];
            else
                disp('\n ERR: proper trajecory not set \n')
            end
        end
    end
end

%% Setting stiffness and damping
if(strcmp(param,'varying') == 1)
    k1 = [repmat(model.k(1),length(tsteady),1),repmat(model.k(2),length(tsteady),1),repmat(model.k(3),length(tsteady),1)];
    k2 = [repmat(0.0,length(tb),1),repmat(0.0,length(tb),1),repmat(0.0,length(tb),1)];
    model.k = [k1; k2];
    
    c1 = [repmat(model.c(1),length(tsteady),1),repmat(model.c(2),length(tsteady),1),repmat(model.c(3),length(tsteady),1)];
    c2 = [repmat(0.0,length(tb),1),repmat(0.0,length(tb),1),repmat(0.0,length(tb),1)];
    model.c = [c1; c2];
else
    
    k1 = [repmat(model.k(1),length(tc),1),repmat(model.k(2),length(tc),1),repmat(model.k(3),length(tc),1)];
    model.k = k1;
    
    c1 = [repmat(model.c(1),length(tc),1),repmat(model.c(2),length(tc),1),repmat(model.c(3),length(tc),1)];
    model.c = c1;
end

model.w = [model.k model.c];

%% Proper initial orientation and world gravity
model.phi0 = Phi(1,:)';
model.B_R_G = euler2dcm(model.phi0);
model.G_g = model.B_R_G'*model.B0_g;

%% Required external wrenches
[f_B1_t, mu_B1_t, f_B2_t, mu_B2_t, omega_B, ddPos, B_g, w] = rigidBodyInvDyn(Phi, dPhi, dPos, ddPos, tc, tDes, model, 0.66, 0.5);


%% System
dPos = @(tDes)interp1(tc,dPos,tDes)';
Phi = @(tDes)interp1(tc,Phi,tDes)';

x0= [dPos(tDes(1));omega_B(tDes(1));zeros(12,1);Phi(tDes(1))];
end