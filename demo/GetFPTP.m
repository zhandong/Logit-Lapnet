function [FPR TPR]=GetFPTP(theta,theta_hat,e)
%% function to get FPR and TPR for simulations and applications
%% ROC is obtained thereby
%% Inputs are theta and estimated theta_hat, e: the tiny threshold
%% Outputs are FPR and TPR that calculated
%% theta is the given coefficients in simulations; theta_hat
%% is the estimated coefficients
thea = abs(theta) > 0;    % transform coefficients to binary values
thea_hat = abs(theta_hat) > e; % convert estimated coefficients to binary values
A = sum(~thea.*~thea_hat);  % A: TN
B = sum(~thea.*thea_hat);   % B: FP
C = sum(thea.*~thea_hat);   % C: FN
D = sum(thea.*thea_hat);    % D: TP
FPR = B/(B+A);         % FPR=FP/(FP+TN)
TPR = D/(D+C);          % TPR=TP/(TP+FN)
