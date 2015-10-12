close all;
clc;
clear;

dataBaseFolder = sprintf('./data/withoutCompliancewithoutSkin/');

load(strcat(dataBaseFolder,'errorResultPeriod005.mat'));
orientationError001 = orientationError;
orientationRMSE001 = orientationRMSE;
model001 = model;
t001 = tKalman;

load(strcat(dataBaseFolder,'errorResultPeriod01.mat'));
orientationError01 = orientationError;
orientationRMSE01 = orientationRMSE;
model01 = model;
t01 = tKalman;

load(strcat(dataBaseFolder,'errorResultPeriod05.mat'));
orientationError1 = orientationError;
orientationRMSE1 = orientationRMSE;
model1 = model;
t1 = tKalman;

tMin = max([t001(1,1);t01(1,1);t1(1,1)]);
tMax = min([t001(end);t01(end);t1(end)]);
dtKalman = min([model001.dtKalman;model01.dtKalman;model1.dtKalman]);
t = linspace(tMin,tMax,(tMax - tMin)/dtKalman);

orientationError001 = interp1(t001,orientationError001,t);
orientationError01 = interp1(t01,orientationError01,t);
orientationError1 = interp1(t1,orientationError1,t);



comp = ['x','y','z'];

%% Orientation RMSE
figure();
stateVar = 19:21;
numSubFig=length(stateVar);
name = {'Orientation RMS Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
    j = stateVar(i);
        subplot(subFig(1),subFig(2),i);
y = [orientationRMSE001(i) ; orientationRMSE01(i); orientationRMSE1(i)];
bar(y);
        xlabel('Period 0.05s       Period 0.1s       Period 0.5s');
        ylabel(strcat('\phi_',comp(i),' m/s'));
        title(name{i});
end
% h = legend('ekf','jekf','dekf');
% h.Box = 'off';
% h.Location = 'northwest';
% h.FontSize = 7;



%% Orientation Error
figure();
stateVar = 19:21;
numSubFig=length(stateVar);
name = {'Orientation Error','',''};
  
subFig = [numSubFig,1];

for i = 1:numSubFig
        
      subplot(subFig(1),subFig(2),i);
        plot(orientationError001(:,i));
        hold on
        plot(orientationError01(:,i));
        hold on
        plot(orientationError1(:,i));
        hold on

        xlabel('time');
        ylabel(strcat('\phi_',comp(i),' m/s'));
        h = legend('t = 0.05s','t = 0.1s','t = 0.5s');

        title(name{i});
end
% h = legend('t = 0.01s','t = 0.1s','t = 1s');
h.Box = 'off';
h.Location = 'northeast';
h.FontSize = 7;

