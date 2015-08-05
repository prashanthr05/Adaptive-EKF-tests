function [kalmanQParams,kalmanRParams,k,c] = setupCovariancesForExperiments();

%kalmanQParams = cell(3);
%% params for experiment 1 
% process model variances
% param ordering : [a_Q, omega_Q,	f_Q,	mu_Q,	phi_Q,	k_Q]
kalmanQParams = [4.0,	10.0,	5.0,	8.0,	0.5,	0.0075];  

% measurement model variances
% param ordering : [f_R,    mu_R,   a_R,    omega_R,    skin_R]
kalmanRParams = [1.5,    2.75,   1.25,   4.5,        25.75];    

% initial stiffness (defined but not used for all experiments)
%k = [k_xx; k_yy; k_zz]
k = [0; 0; 0];

% initial damping (defined but not used for all experiments)
%c = [c_xx; c_yy; c_zz]
c = [0; 0; 0];


end
