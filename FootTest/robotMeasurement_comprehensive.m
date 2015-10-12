function [yMeas,tMeas,model] = robotMeasurement_comprehensive(dtKalman, model, plots, t_min, t_max, measurementType,processType,numberOfExperiment ,whichLeg, whichSkin,dataSet,expt)

if(nargin<1)
    measurementType = 'withoutSkin';
    dtKalman = 0.01;
    model = struct();
    plots = 0;
    t_min = 0;
    t_max = inf;
    dataSet = 'BWD'
    numberOfExperiment = 1;
    expt = 'leg';
end

if (nargin<6)
    measurementType = 'withoutSkin';
    numberOfExperiment = 1;
    whichLeg  = 'right';
    whichSkin = 'right';
    dataSet = 'BWD';
end

%% Setting up fake measurements for check
fakeIMU = 'false';
fakeFT_f = 'false';
fakeFT_mu = 'false';

fakeAccl = [0;0;-9.8];
fakeMu_leg = [0;0;0];
fakeF_leg = [0;0;0];

[model,tMax,leg_ft,foot_ft,inertial,skin,transforms] = robotSensorParams_comprehensive(whichLeg,whichSkin,numberOfExperiment,t_min,t_max,dtKalman,model,dataSet,expt,measurementType);

model.dtKalman = dtKalman;
t = linspace(t_min,tMax,(tMax - t_min)/dtKalman);
tCalib = linspace(0,t_min,(t_min)/dtKalman);

%% Calibrating the measurements
a_omega_calib = interp1(inertial.t,inertial.data,tCalib);

omega_calib = a_omega_calib(:,7:9);
omegaOffset = mean((omega_calib'),2);
phi_pre_calib = (transforms.B_R_imu*a_omega_calib(:,1:3)');
phiMean = mean(phi_pre_calib,2);


%% raw interpolated data
f_leg_raw = interp1(leg_ft.t,leg_ft.f,t);
mu_leg_raw = interp1(leg_ft.t,leg_ft.mu,t);
f_foot_raw = interp1(foot_ft.t,foot_ft.f,t);
mu_foot_raw = interp1(foot_ft.t,foot_ft.mu,t);

a_omega_raw = interp1(inertial.t,inertial.data,t);
a_raw = a_omega_raw(:,4:6);
omega_raw = a_omega_raw(:,7:9);

phi = (transforms.B_R_imu*a_omega_raw(:,1:3)');
phi = phi - repmat(phiMean,1,size(phi,2));
model.phi_groundtruth = repmat(model.phi0.*(180/pi),1,size(a_raw,1)) - phi;

delta_raw = interp1(skin.t,skin.data,t);

rawData.f_leg = f_leg_raw;
rawData.mu_leg = mu_leg_raw;
rawData.f_foot = f_foot_raw;
rawData.mu_foot = mu_foot_raw;
rawData.a = a_raw;
rawData.omega = omega_raw;

% Total normal force from the skin
delta = dataPostProcessing(delta_raw,'normalForces');
fc_z = computeTotalForce(delta,'normalForces');


%IMU
a = transforms.B_R_imu*a_raw';
omegaCentered = (omega_raw')' - repmat(omegaOffset',size(omega_raw,1),1);
omega = (transforms.B_R_imu*omegaCentered')'.*(pi/180);


        f_mu_leg = transforms.B_adjT_ankle*[f_leg_raw';mu_leg_raw'];
        f_mu_foot = transforms.B_adjT_ankle*[f_foot_raw';mu_foot_raw'];
        
        f_leg = f_mu_leg(1:3,:) ;%- f_calib_delta*ones(1,length(t));
        mu_leg = f_mu_leg(4:6,:) ;%- mu_calib_delta*ones(1,length(t));
        f_foot = f_mu_foot(1:3,:) ;%- f_calib_delta*ones(1,length(t));
        mu_foot = f_mu_foot(4:6,:) ;%- mu_calib_delta*ones(1,length(t));
   
%% Fake Measurements   
   if(strcmp(fakeIMU,'true')==1)
   	aPerfect = fakeAccl*ones(1,size(a,2));
   	omegaPerfect = zeros(size(omega));
   	a = aPerfect;
   	omega = omegaPerfect;
   end
   
   
   if(strcmp(fakeFT_f,'true')==1)
   	flegPerfect = fakeF_leg*ones(1,size(f_leg,2));
   	fgPerfect = model.m*transpose(euler2dcm(model.phi0))*model.G_g*ones(1,size(f_foot,2));
   	ffootPerfect = - 0.5*fgPerfect;
   	flegPerfect = -ffootPerfect;
   	fczPerfect = ffootPerfect(3,:);
   	f_leg = flegPerfect;
   	f_foot = ffootPerfect;
   	fc_z = fczPerfect;
   end
   
      if(strcmp(fakeFT_mu,'true')==1)
      	mulegPerfect = fakeMu_leg*ones(1,size(mu_leg,2));
      	mufootPerfect = fakeMu_leg*ones(1,size(mu_foot,2));      
	
	mu_leg = mulegPerfect;
	mu_foot = mufootPerfect;     
     end
        
%%   Seting up measurement and initial conditions

fo = f_foot;
muo = mu_foot;
fc = fo + model.m*a;
muc = muo;



        fprintf('\n setting up robot measurements for leg\n');
        if(strcmp(processType,'withoutCompliance')==1)
            model.x0 = [zeros(3,1);zeros(3,1);fo(:,1);muo(:,1);fc(:,1);muc(:,1);model.phi0];
        else if(strcmp(processType,'withCompliance')==1)
            model.x0 = [zeros(3,1);zeros(3,1);fo(:,1);muo(:,1);fc(:,1);muc(:,1);model.phi0;model.k;model.c];
            end
        end
    
    	if(strcmp(measurementType,'withoutSkin')==1)
	        yMeas = [a;omega';fo;muo;fc;muc]';
        else
        	yMeas = [a;omega';fo;muo;fc;muc;fc_z']';
        end

tMeas = t;
if (strcmp(model.measurementPlots,'makePlots')==1)
    figure(6);
    %
    subplot(2,2,1);
    plot(t,a_raw);
    xlabel('time (sec)');
    ylabel('a - g (m/s^2)');
    legend('aX', 'aY', 'aZ');
    axis tight;
    title('Raw Proper Acceleration');
    
    subplot(2,2,2);
    plot(t,omega_raw);
    xlabel('time (sec)');
    ylabel('\omega (deg/s)');
    legend('\omega_{X}', '\omega_{Y}', '\omega_{Z}');
    axis tight;
    title('Raw Angular Velocity');
    
    subplot(2,2,3);
    plot(t,a);
    xlabel('time (sec)');
    ylabel('a - g (m/s^2)');
    legend('aX', 'aY', 'aZ');
    axis tight;
    title('Proper Acceleration expressed in Body frame');
    
    subplot(2,2,4);
    plot(t,omega);
    xlabel('time (sec)');
    ylabel('\omega (deg/s)');
    legend('\omega_{X}', '\omega_{Y}', '\omega_{Z}');
    axis tight;
    title('Angular Velocity expressed in Body frame');
    
    
    figure(7);
    %
    subplot(2,2,1);
    plot(t,fo);
    xlabel('time (sec)');
    ylabel('fo (N)');
    legend('foX', 'foY', 'foZ');
    axis tight;
    title('Body Force expressed in Body Frame');
    
    subplot(2,2,2);
    plot(t,muo);
    xlabel('time (sec)');
    ylabel('\mu_{o} (Nm)');
    legend('\muo_{X}', '\muo_{Y}', '\muo_{Z}');
    axis tight;
    title('Body Moment in Body Frame');
    
    subplot(2,2,3);
    plot(t,fc);
    xlabel('time (sec)');
    ylabel('fc (Nm)');
    legend('fcX', 'fcY', 'fcZ');
    axis tight;
    title('Contact Force expressed in Body frame');
    
    subplot(2,2,4);
    plot(t,muc);
    xlabel('time (sec)');
    ylabel('\mu_{c} (Nm)');
    legend('\muc_{X}', '\muc_{Y}', '\muc_{Z}');
    axis tight;
    title('Contact expressed in Body frame');
    
end

end
