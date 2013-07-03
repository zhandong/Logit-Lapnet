%% script to simulate testing and training so as to get the squared mean of
%% the error, i.e., the PMSE (Penalized Mean Squared Errors), scripts 
%% contain three parts that for Elastic, Lasso and LogitNet respectively:
%% getPMSEENet, getPMSELasso, getPMSE
%% the last one, say, this script is to get PMSE for logistic-lapnet 
clear all;
load coeff1;  %  coefficients in simulation1
% load coeff2; % coefficients in simulation2
% load coeff3; % coefficients in simulation3
% load coeff4; % coefficients in simulation4
lambda_opt=.00162975083; % optimal lambda that obtained by CV
alpha1=.2;
fprintf('\n');
disp('--- Cross Validation is in process, interface is kinda wordy though ..... \n');
fprintf('\n');
for i=1:50
    [Xm,y,L]=SimNetData(beta1,200,100);
    theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1));
    [Xm1,y1,L1]=SimNetData(beta1,200,100);
    logisError(i)=logisSimError(Xm1,y1,theta_hat);
end
a=mse(logisError);
a1=mean(logisError);
fprintf('\n The mean squared error is: %15.6f \n',a)
standerr=sqrt(sum((logisError-a1).^2));
fprintf('The standard error is: %15.6f \n',standerr)
save PMSEAll % save work space