function [kalmanQParams,kalmanRParams,k,c] = gazeboCovariancesForExperiments(filter);

%kalmanQParams = cell(3);
%% params for experiment 1 
% process model variances
% param ordering : [a_Q, omega_Q,	f_Q,	mu_Q,	phi_Q,	k_Q, c_Q]
% kalmanQParams = [4.0,	10.0,	5.0,	8.0,	0.5,	0.025, 0.0025];  
% kalmanQParams = [4.0,	10.0,	5.0,	2.0,	0.5,	0.025, 0.0025];  

% kalmanQParams = [0.5,	0.5,	2.0,	2.0,	0.5,	0.2595, 0.025];

% kalmanQParams = [0.5,	0.5,	2.0,	2.0,	0.5,	0.002595, 0.00025];

% kalmanQParams = [0.5,	0.5,	0.5,	0.5,	1.625,	0.02595, 0.0025];
% kalmanQParams = [0.0,	0.5,	1.5,	0.5,	0.0,	0.0, 0.00];


kalmanQParams = [4.0,	2.0,	5.0,	0.5,	0.0,	0.025, 0.0025];  


% kalmanQParams = [0.0,	0.005,	0.5,	0.5,	0.0025,	0.00002595, 0.00025];
% kalmanQParams = [0.00005,	0.00005,	0.0005,	0.0005,	0.0005,	0.0005, 0.0005];



% measurement model variances
% param ordering : [f_R,    mu_R,   a_R,    omega_R,    skin_R]
kalmanRParams = [1.5,    2.75,   1.25,   4.5,        25.75];
% kalmanRParams = [3.0901,    0.0188,   1.25,   0.3080,        25.75];

% kalmanRParams = [1.5,    2.75,   1.25,   0.5,        25.75];

% kalmanRParams = [0.5,    0.5,   0.5,   0.5,        0.5];
% kalmanRParams = [0.6348,    0.0278,   6.5287e-04,   1.0003e-04,        25.75];

% kalmanRParams = [2.5463,    0.0178,   1.25,   0.0174,        25.75];

% initial stiffness (defined but not used for all experiments)
%k = [k_xx; k_yy; k_zz]
k = [0.05; 0.05; 0.05];
% k = [0.25; 0.25; 0.25];
% initial damping (defined but not used for all experiments)
%c = [c_xx; c_yy; c_zz]
c = [0.005; 0.005; 0.005];


end
