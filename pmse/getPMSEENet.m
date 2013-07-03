%% This script is to simulate the 3rd test for Elastic net 50 times and get the error mean
clear all;
load coeff1; % different coefficients in four simulations
% load coeff2;
% load coeff3;
% load coeff4;
lambda_opt=49.77023564; %% the optimal lambda for ENet got by CV
alpha1=.2;
logisError=zeros(50,1);
Sp=zeros(50,1);
Sen=zeros(50,1);
disp('--- Data is in process, wait patiently and enjoy the wordy interface ....');
fprint('\n')
for i=1:50
    [Xm,y,L]=SimEnetData(beta1,200,100); %% get simulation training data for ENet
    theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1)); 
    [Sp(i) Sen(i)]=GetSpSen(beta1,theta_hat, 0.0001); %% obtain specificity and sensitivity
    [Xm1,y1,L1]=SimEnetData(beta1,200,100); %% testing data
    logisError(i)=logisSimError(Xm1,y1,theta_hat); %% store the error for each time
end

%% PMSE and the related standard error
squaremean=mse(logisError);
a1=mean(logisError);
fprintf('\nThe PMSE is: %15.6f \n',squaremean)
standerr=sqrt(sum((logisError-a1).^2))/sqrt(50);
fprintf('The standard error of PMSE is: %15.6f \n',standerr)
%% sensitivity and the corresponding standard error
a1=mean(Sen);
fprintf('The sensitivity is: %15.6f \n',a1)
standerrsen=sqrt(sum((Sen-a1).^2))/sqrt(50);
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)
%% specificity and related standard error
specificity=mean(Sp);
fprintf('The spacificity is: %15.6f \n',specificity)
standerrsp=sqrt(sum((Sp-specificity).^2))/sqrt(50);
fprintf('The standard error is of spacificity: %15.6f \n',standerrsp)
save PMSEENet3 % save work space
