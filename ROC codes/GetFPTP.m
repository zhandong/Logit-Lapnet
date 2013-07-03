function [FPR TPR]=GetFPTP(theta,theta_hat,e)
%% function to get FPR (False Positive Rate) and TPR (True Positive Rate) 
%% for simulations and applications
%% ROC is obtained thereby
%% Inputs are objective coefficients theta and estimated coefficients 
%% theta_hat, e: tiny threshold
%% Outputs are FPR and TPR that calculated
thea = abs(theta) > 0;    % transform coefficients to binary values
thea_hat = abs(theta_hat) > e; % convert estimated coefficients to binary values
A = sum(~thea.*~thea_hat);  % A: TN
B = sum(~thea.*thea_hat);   % B: FP
C = sum(thea.*~thea_hat);   % C: FN
D = sum(thea.*thea_hat);    % D: TP
FPR = B/(B+A);         % FPR=FP/(FP+TN)
TPR = D/(D+C);          % TPR=TP/(TP+FN)
