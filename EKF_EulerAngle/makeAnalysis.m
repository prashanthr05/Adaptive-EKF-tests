close all;
clc;
clear;

idx = 1;

dataBaseFolder = sprintf('./data/');

load(strcat(dataBaseFolder,'analysis.mat'));
tK = tKalman';
tc = 800;
comp = ['x','y','z'];
numPlots = size(stateError,1);
%% Orientation error
figure();

stateVar = 19:21;
numSubFig=length(stateVar);
name = {'Orientation Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        subplot(subFig(1),subFig(2),i);
        
     for j = idx : numPlots  
         theta = stateError{j}(1:tc,stateVar);
%          theta = XUpt{j}(:,stateVar);
        plot(tK(:,idx:tc),theta(idx:end,i));
        xlabel('Time (s)');
        ylabel(strcat('\phi_',comp(i),' degs'));
        hold on    
     end    
    title(name{i});
    
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'southwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/analysis/orientationError','-opengl');


%% Linear Velocity error
figure();

stateVar = 1:3;

numSubFig=length(stateVar);
name = {'Linear Velocity Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        subplot(subFig(1),subFig(2),i);
        
     for j = idx : numPlots  
         v = stateError{j}(1:tc,stateVar);
        plot(tK(:,idx:tc),v(idx:end,i));
        xlabel('Time (s)');
        ylabel(strcat('v_',comp(i),' m/s'));
        hold on    
     end    
    title(name{i});
    
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'southwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/analysis/linVelError','-opengl');

%% Angular Velocity error
figure();

stateVar = 4:6;
numSubFig=length(stateVar);
name = {'Angular Velocity Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        subplot(subFig(1),subFig(2),i);
        
     for j = idx : numPlots  
         omega = stateError{j}(1:tc,stateVar);
        plot(tK(:,idx:tc),omega(idx:end,i));
        xlabel('Time (s)');
        ylabel(strcat('\omega_',comp(i),' rad/s'));
        
        hold on    
     end    
    title(name{i});
    
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'southwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/analysis/angularVelError','-opengl');



%% Stiffness error
figure();
paramVar = 1:3;
numSubFig=length(paramVar);
name = {'Stiffness Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        subplot(subFig(1),subFig(2),i);
        
     for j = idx : numPlots  
         K = paramError{j}(1:tc,paramVar);
%          theta = XUpt{j}(:,stateVar);
        plot(tK(:,idx:tc),K(idx:end,i));
        xlabel('Time (s)');
        ylabel(strcat('K_',comp(i),' Nm/rad'));
        hold on    
     end    
    title(name{i});
    
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'southwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/analysis/stiffnessError','-opengl');



%% Damping error

figure();

paramVar = 4:6;
numSubFig=length(paramVar);
name = {'Damping Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        subplot(subFig(1),subFig(2),i);
        
     for j = idx : numPlots  
         C = paramError{j}(1:tc,paramVar);
%          theta = XUpt{j}(:,stateVar);
        plot(tK(:,idx:tc),C(idx:end,i));
        xlabel('Time (s)');
        ylabel(strcat('C_',comp(i),' Nms/rad'));
        hold on    
     end    
    title(name{i});
    
end
h = legend('ekf','jekf','dekf');
h.Box = 'off';
h.Location = 'southwest';
h.FontSize = 7;
set(gca,'FontSize',12);
set(gcf,'Renderer','OpenGL');
print('-djpeg','-r200','./plots/analysis/dampingError','-opengl');




