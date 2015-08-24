clc
clear;
close all;

%% Add paths
utilities    = genpath('./utils');
dynFunctions = genpath('./dynamicsFunctions');
addpath(utilities,dynFunctions);

%% Load data
dataFolder = './Measurement/';
load(strcat(dataFolder,'measurement.mat'));
dtForDyn = 1e-3;
tc = model.t_min : dtForDyn : model.t_max;

model.invDynPlots = 'makePlots';

fo = sysX(7:9,:);
muo = sysX(10:12,:);
fc = sysX(13:15,:);
muc = sysX(16:18,:);

f_B_o = @(tc)interp1(tDes,fo',tc)';
f_B_c = @(tc)interp1(tDes,fc',tc)';
mu_B_o = @(tc)interp1(tDes,muo',tc)';
mu_B_c = @(tc)interp1(tDes,muc',tc)';


sysX = interp1(tDes,sysX',tc);
p = model;

m = p.m;
I_B = p.I;
phi0 = p.phi0;
G_g = p.G_g;
k = p.k;
K = diag([k(1),k(2),k(3)]);
c = p.c;
C = diag([c(1),c(2),c(3)]);


z0 = p.x0;

x0 = [z0(1:3);z0(4:6);z0(19:21)]';
tSpan = tc;


% rigid = @(t,x)[-S(x(4:6)) * x(1:3) + 1/m * f_B_o(t) - 1/m*f_B_c(t) + euler2dcm(x(7:9))*G_g;...
%                I_B \ (-S(x(4:6)) * (I_B * x(4:6)) + mu_B_o(t) - mu_B_c(t) - (- K'*K*(x(7:9) - phi0) - C*x(4:6))); ...
%                (cell2mat(Tomega_dphi(x(7:9)','t')))\x(4:6)];

rigid = @(t,x)integrateForward(t,x,f_B_o(t),mu_B_o(t),f_B_c(t),mu_B_c(t),p);
% odeSettings = odeset('InitialStep', 1e-8, 'MaxStep',1e-4);

[t,x] = ode45(rigid,tSpan,x0);%,odeSettings);

if(strcmp(model.invDynPlots,'makePlots') == 1)
    figure(1);
    subplot(3,1,1);
    plot(tc,x(:,1)); hold on;
    plot(tc,sysX(:,1)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,1);sysX(:,1)])-1 max([x(:,1);sysX(:,1)])+1]);
    legend('forDyn','invDyn');
    title('Linear velociity vx')
    subplot(3,1,2);
    plot(tc,x(:,2)); hold on;
    plot(tc,sysX(:,2)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,2);sysX(:,2)])-1 max([x(:,2);sysX(:,2)])+1]);
    legend('forDyn','invDyn');
    title('Linear velociity vy')
    subplot(3,1,3);
    plot(tc,x(:,3)); hold on;
    plot(tc,sysX(:,3)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,3);sysX(:,3)])-1 max([x(:,3);sysX(:,3)])+1]);
    legend('forDyn','invDyn');
    title('Linear velociity vz');
        set(gca,'FontSize',12);
    set(gcf,'Renderer','OpenGL');
    print('-depsc2','-r200','./plots/InvDyn/LinVelTraj','-opengl');
    
    
    
    figure(2);
    subplot(3,1,1);
    plot(tc,x(:,4)); hold on;
    plot(tc,sysX(:,4)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,4);sysX(:,4)])-1 max([x(:,4);sysX(:,4)])+1]);
    legend('forDyn','invDyn');
    title('Angular velociity vx')
    subplot(3,1,2);
    plot(tc,x(:,5)); hold on;
    plot(tc,sysX(:,5)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,5);sysX(:,5)])-1 max([x(:,5);sysX(:,5)])+1]);
    legend('forDyn','invDyn');
    title('Angular velociity vy')
    subplot(3,1,3);
    plot(tc,x(:,6)); hold on;
    plot(tc,sysX(:,6)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,6);sysX(:,6)])-1 max([x(:,6);sysX(:,6)])+1]);
    legend('forDyn','invDyn');
    title('Angular velociity vz');
    set(gca,'FontSize',12);
    set(gcf,'Renderer','OpenGL');
    print('-depsc2','-r200','./plots/InvDyn/AngVelTraj','-opengl');
    
    figure(3);
    subplot(3,1,1);
    plot(tc,x(:,7)); hold on;
    plot(tc,sysX(:,19)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,7);sysX(:,19)])-1 max([x(:,7);sysX(:,19)])+1]);
    legend('forDyn','invDyn');
    title('\phi_{x}')
    subplot(3,1,2);
    plot(tc,x(:,8)); hold on;
    plot(tc,sysX(:,20)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,8);sysX(:,20)])-1 max([x(:,8);sysX(:,20)])+1]);
    legend('forDyn','invDyn');
    title('\phi_{y}')
    subplot(3,1,3);
    plot(tc,x(:,9)); hold on;
    plot(tc,sysX(:,21)','--m'); hold on;
    axis([p.t_min p.t_max min([x(:,9);sysX(:,21)])-1 max([x(:,9);sysX(:,21)])+1]);
    legend('forDyn','invDyn');
    title('\phi_{z}');
    set(gca,'FontSize',12);
    set(gcf,'Renderer','OpenGL');
    print('-depsc2','-r200','./plots/InvDyn/OrientationTraj','-opengl');
    
    
end
