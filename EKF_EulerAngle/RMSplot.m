close all;
clc;
clear;

dataBaseFolder = sprintf('./data/RMSE/');

load(strcat(dataBaseFolder,'constantParam.mat'));
stateRMSEconstant = stateRMSE;
paramRMSEconstant = paramRRSE;

load(strcat(dataBaseFolder,'singlevaryingParam.mat'));
stateRMSE1D = stateRMSE;
paramRMSE1D = paramRRSE;

load(strcat(dataBaseFolder,'multivaryingParam.mat'));
stateRMSE3D = stateRMSE;
paramRMSE3D = paramRRSE;

comp = ['x','y','z'];

%% Stiffness RMSE
figure();
paramVar = 1:3;
numSubFig=length(paramVar);
name = {'Stiffness RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
    j = paramVar(i);
        subplot(subFig(1),subFig(2),i);
y = [paramRMSEconstant{2}(j) paramRMSEconstant{3}(j); paramRMSE1D{2}(j) paramRMSE1D{3}(j); paramRMSE3D{2}(j) paramRMSE3D{3}(j)];
bar(y);
        xlabel('Constant Params       1-D Varying Params       3-D varying Params');
        ylabel(strcat('K_',comp(i),' Nm/rad'));
        title(name{i});
end
h = legend('jekf','dekf');
h.Box = 'off';
h.Location = 'northwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/performance/stiffnessRMSError','-opengl');

%% Damping RMSE
figure();
paramVar = 4:6;
numSubFig=length(paramVar);
name = {'Damping RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
    j = paramVar(i);
        subplot(subFig(1),subFig(2),i);
y = [paramRMSEconstant{2}(j) paramRMSEconstant{3}(j); paramRMSE1D{2}(j) paramRMSE1D{3}(j); paramRMSE3D{2}(j) paramRMSE3D{3}(j)];
bar(y);
        xlabel('Constant Params       1-D Varying Params       3-D varying Params');
        ylabel(strcat('C_',comp(i),' Nms/rad'));
        title(name{i});
end
h = legend('jekf','dekf');
h.Box = 'off';
h.Location = 'northwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/performance/dampingRMSError','-opengl');



%% Lin Vel RMSE
figure();
stateVar = 1:3;
numSubFig=length(stateVar);
name = {'Linear Velocity RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
    j = stateVar(i);
        subplot(subFig(1),subFig(2),i);
y = [stateRMSEconstant{1}(j) stateRMSEconstant{2}(j) stateRMSEconstant{3}(j); stateRMSE1D{1}(j) stateRMSE1D{2}(j) stateRMSE1D{3}(j); stateRMSE3D{1}(j) stateRMSE3D{2}(j) stateRMSE3D{3}(j)];
bar(y);
        xlabel('Constant Params       1-D Varying Params       3-D varying Params');
        ylabel(strcat('v_',comp(i),' m/s'));
        title(name{i});
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'northwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/performance/linvelRMSError','-opengl');

%% Ang Vel RMSE
figure();
stateVar = 4:6;
numSubFig=length(stateVar);
name = {'Angular Velocity RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
    j = stateVar(i);
        subplot(subFig(1),subFig(2),i);
y = [stateRMSEconstant{1}(j) stateRMSEconstant{2}(j) stateRMSEconstant{3}(j); stateRMSE1D{1}(j) stateRMSE1D{2}(j) stateRMSE1D{3}(j); stateRMSE3D{1}(j) stateRMSE3D{2}(j) stateRMSE3D{3}(j)];
bar(y);
        xlabel('Constant Params       1-D Varying Params       3-D varying Params');
        ylabel(strcat('\omega_',comp(i),' rad/s'));
        title(name{i});
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'northwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/performance/angvelRMSError','-opengl');

%% Orientation RMSE
figure();
stateVar = 19:21;
numSubFig=length(stateVar);
name = {'Orientation RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
    j = stateVar(i);
        subplot(subFig(1),subFig(2),i);
y = [stateRMSEconstant{1}(j) stateRMSEconstant{2}(j) stateRMSEconstant{3}(j); stateRMSE1D{1}(j) stateRMSE1D{2}(j) stateRMSE1D{3}(j); stateRMSE3D{1}(j) stateRMSE3D{2}(j) stateRMSE3D{3}(j)];
bar(y);
        xlabel('Constant Params       1-D Varying Params       3-D varying Params');
        ylabel(strcat('\phi_',comp(i),' m/s'));
        title(name{i});
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'northwest';
h.FontSize = 7;

set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/performance/orientationRMSError','-opengl');