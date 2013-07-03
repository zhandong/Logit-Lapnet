%% script to simulate experiments and get the squared mean error 
clc;
clear all;
load coeff1; % coefficients in simulation1
% load coeff2;
% load coeff3;
% load coeff4;
lambda_opt=37.64935807; % optimal lambda for Logit-LapNet, got via CV
alpha1=.2;
disp('--- Data is in process, enjoy the wordy interface  .....');
fprint('\n')
%% vectors to store error and sensitivity, specificity
logisError=zeros(50,1);
Sp=zeros(50,1); % array to store specificity
Sen=zeros(50,1); % array to store sensitivity
for i=1:50
    [Xm,y,L]=SimNetData(beta1,200,100); %% get training set
    theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1));  % estimated theta_hat
    [Sp(i) Sen(i)]=GetSpSen(beta1,theta_hat, 0.0001); % sensitivity and specificity
    [Xm1,y1,L1]=SimNetData(beta1,200,100); %% testing data set 
    logisError(i)=logisSimError(Xm1,y1,theta_hat); % the PMSE
end
fprint('\n')
%% Print the results 
squaremean=mse(logisError);
a1=mean(logisError);
fprintf('\n The PMSE is: %15.6f \n',squaremean)
standerr=sqrt(sum((logisError-a1).^2))/sqrt(50);
fprintf('The standard error of PMSE is: %15.6f \n',standerr)
%% sensitivity and standard error
sensitivity=mean(Sen);
fprintf('The sensitivity is: %15.6f \n',sensitivity)
standerrsen=sqrt(sum((Sen-sensitivity).^2))/sqrt(50);
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)
%% specificity and standard error
specificity=mean(Sp);
fprintf('The specificity is: %15.6f \n',specificity)
standerrsp=sqrt(sum((Sp-specificity).^2))/sqrt(50);
fprintf('The standard error of spacificity: %15.6f \n',standerrsp)
%% Save the results
save PMSENet 

