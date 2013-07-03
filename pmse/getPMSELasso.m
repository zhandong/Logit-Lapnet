%% This script is to simulate the first Lasso test 50 times and get the error mean
clear all;
load coeff1;  % different coefficients in four simulations
% load coeff2; % coefficients in simulation2
% load coeff3; % coefficients in simulation3
% load coeff4; % coefficients in simulation4
lambda_opt=.00162975083; % optimal lambda got by CV
alpha1=1;
logisError=zeros(50,1);
Sp=zeros(50,1);
Sen=zeros(50,1);
disp('--- Data is in process, wait patiently and enjoy the wordy interface ....');
fprint('\n')
for i=1:50
    [Xm,y,L]=SimLassoData(beta1,200,100); %% get the training data for Lasso simulation
    theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1));
    [Sp(i) Sen(i)]=GetSpSen(beta1,theta_hat, 0.0001); %% sensitivity and specificity
    [Xm1,y1,L1]=SimNetData(beta1,200,100); %% testing data
    logisError(i)=logisSimError(Xm1,y1,theta_hat); %% errors
end

%% PMSE and standard errors 
squaremean=mse(logisError);
a1=mean(logisError);
fprintf('The PMSE is: %15.6f \n',squaremean)
standerr=sqrt(sum((logisError-a1).^2))/sqrt(50);
fprintf('The standard error of PMSE is: %15.6f \n',standerr)

%% sensitivity and standard error for it
sensitivity=mean(Sen);
fprintf('The sensitivity is: %15.6f \n',sensitivity)
standerrsen=sqrt(sum((Sen-sensitivity).^2))/sqrt(50);
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)

%% specificity and standard error of it
specificity=mean(Sp);
fprintf('The specificity is: %15.6f \n',specificity)
standerrsp=sqrt(sum((Sp-specificity).^2))/sqrt(50);
fprintf('The standard error is of spacificity: %15.6f \n',standerrsp)

save PMSELasso 
