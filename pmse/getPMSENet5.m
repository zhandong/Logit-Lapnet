%% This script is to simulate the fourth model 50 times and get the mean of the error 
%% It is for Net and get the related PMSE squared mean and standard error
clear all;
beta1=zeros(2200,1);  % beta is a matlab function , so to use beta1 here
beta1(1)=5;
beta1(12)=-5;
beta1(23)=3;
beta1(34)=-3;
beta1(5:11)=3/(sqrt(10))*zeros(7,1);
beta1(2:4)=-5/(sqrt(10))*ones(3,1);
beta1(16:22)=-3/(sqrt(10))*zeros(7,1);   % there are more zero coefficients in fifth model
beta1(13:15)=5/(sqrt(10))*ones(3,1);
beta1(24:26)=-3/(sqrt(10))*ones(3,1);
beta1(27:33)=5/(sqrt(10))*zeros(7,1);
beta1(35:37)=3/(sqrt(10))*ones(3,1); 
beta1(38:44)=-5/(sqrt(10))*zeros(7,1);  % to generate the beta in simulation

lambda_opt=65.793322; % the optimal lambda for Net, which is obtained via CV 
alpha1=.2;   % alpha is set to .2 here.
logisError=zeros(50,1);
Sp=logisError;
Sen=logisError;

for i=1:50
    [Xm,y,L]=SimNetData(beta1,200,100); %% get the training data 
    theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1)); %% the theta_hat obtained via the logistic regression model
    [Sp(i) Sen(i)]=GetSpSen(beta1,theta_hat, 0.0001); % Get sensitivity and specificity for each training data
    [Xm1,y1,L1]=SimNetData(beta1,200,100); %% another different set of data used as testing data 
    logisError(i)=logisSimError(Xm1,y1,theta_hat); % Get the error of y and y_hat
end


squaremean=mse(logisError); % squared mean of the errors
a1=mean(logisError); % mean of the errors
fprintf('The PMSE is: %15.6f \n',squaremean)% Print the results out
standerr=sqrt(sum((logisError-a1).^2))/sqrt(50); % standard errors
fprintf('The standard error of PMSE is: %15.6f \n',standerr)


sensitivity=mean(Sen); % take the mean of sensitivities as the sensitivity
fprintf('The sensitivity is: %15.6f \n',sensitivity)
standerrsen=sqrt(sum((Sen-sensitivity).^2))/sqrt(50); % the standard error of the sensitivity
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)


spacificity=mean(Sp); % the specificity that calculated
fprintf('The spacificity is: %15.6f \n',spacificity)
standerrsp=sqrt(sum((Sp-spacificity).^2))/sqrt(50);  % standard error of the specificity
fprintf('The standard error is of spacificity: %15.6f \n',standerrsp)

%% Save the results
save PMSENet5 


