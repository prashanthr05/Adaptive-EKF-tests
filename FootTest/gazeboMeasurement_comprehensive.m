function [yMeas,tMeas,model] = gazeboMeasurement_comprehensive(dtKalman, model, plots, t_min, t_max, measurementType,processType,numberOfExperiment ,whichLeg, whichSkin,dataSet,expt)

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

% fakeFlags.fakeIMU = 'false';
% fakeFlags.fakeFT_f = 'false';
% fakeFlags.fakeFT_mu = 'false';
%
% fakeValues.fakeAccl =[ 0 ; 0 ; 9.8];
% fakeValues.fakeF_o = [0;0;0];
% fakeValues.fakeMu_o = [0;0;0];

[model,tMax,leg_ft,foot_ft,inertial,transforms] = gazeboSensorParams_comprehensive(whichLeg,whichSkin,numberOfExperiment,t_max,dtKalman,model,dataSet,expt);

model.dtKalman = dtKalman;
t = linspace(t_min,tMax,(tMax - t_min)/dtKalman);
tCalib = linspace(0,t_min,(t_min)/dtKalman);

a_omega_pre_calib = interp1(inertial.t,inertial.data,tCalib);

f_foot_pre_calib = interp1(foot_ft.t,foot_ft.f,tCalib);
mu_foot_pre_calib = interp1(foot_ft.t,foot_ft.mu,tCalib);

f_mu_foot_pre_calib = transforms.B_adjT_ankle*[f_foot_pre_calib';mu_foot_pre_calib'];

f_foot_pre_calib = f_mu_foot_pre_calib(1:3,:);
mu_foot_pre_calib = f_mu_foot_pre_calib(4:6,:);

a_pre_calib = a_omega_pre_calib(:,4:6);
omega_pre_calib = a_omega_pre_calib(:,7:9);

phi_pre_calib = (transforms.B_R_imu*a_omega_pre_calib(:,1:3)');

fMean = mean(f_foot_pre_calib,2);
muMean = mean(mu_foot_pre_calib,2);
aMean = mean((a_pre_calib'),2);
omegaMean = mean((omega_pre_calib'),2);

phiMean = mean(phi_pre_calib,2);

aDeviation = a_pre_calib' - repmat(aMean,1,size(a_pre_calib,1));
omegaDeviation = omega_pre_calib' - repmat(omegaMean,1,size(omega_pre_calib,1));
fDeviation = f_foot_pre_calib - repmat(fMean,1,size(f_foot_pre_calib,2));
muDeviation = mu_foot_pre_calib - repmat(muMean,1,size(mu_foot_pre_calib,2));
a_R = 0;
omega_R = 0;
f_R = 0;
mu_R = 0;
for i = 1 : length(tCalib)
    a_R = a_R + transpose(aDeviation(1:3,i))*aDeviation(1:3,i);
    omega_R = omega_R + transpose(omegaDeviation(1:3,i))*omegaDeviation(1:3,i);
    f_R = f_R + transpose(fDeviation(1:3,i))*fDeviation(1:3,i);
    mu_R = mu_R + transpose(muDeviation(1:3,i))*muDeviation(1:3,i);
end
omega_R = omega_R/length(tCalib);
a_R = a_R/length(tCalib);
f_R = f_R/length(tCalib);
mu_R = mu_R/length(tCalib);



% raw interpolated data
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


a = transforms.B_R_imu*a_raw';
omegaCentered = (omega_raw')' - repmat(omegaMean',size(omega_raw,1),1);
omega = (transforms.B_R_imu*omegaCentered')'.*(pi/180);

f_mu_leg = transforms.B_adjT_ankle*[f_leg_raw';mu_leg_raw'];
f_mu_foot = transforms.B_adjT_ankle*[f_foot_raw';mu_foot_raw'];

f_leg = f_mu_leg(1:3,:);
mu_leg = f_mu_leg(4:6,:);
f_foot = f_mu_foot(1:3,:);
mu_foot = f_mu_foot(4:6,:);

fo = f_foot;
muo = mu_foot;
fc = fo + model.m*a;
muc = muo;

yMeas = [a;omega';fo;muo;fc;muc]';
if(strcmp(processType,'withoutCompliance')==1)
    model.x0 = [zeros(3,1);zeros(3,1);fo(:,1);muo(:,1);fc(:,1);muc(:,1);model.phi0];
else if(strcmp(processType,'withCompliance')==1)
        model.x0 = [zeros(3,1);zeros(3,1);fo(:,1);muo(:,1);fc(:,1);muc(:,1);model.phi0;model.k;model.c];
    end
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
