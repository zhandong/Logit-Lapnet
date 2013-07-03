clear all; clc
load coeff1; % the first simulation coefficients
% load coeff2;  % the 2nd simulation
%  load coeff3; % the 3rd simulation
%  load coeff4; % the 4th simulation
 %% Xm, y, L in simulation
[Xm,y,L]=SimLassoData(beta1,200,100); % simulate lasso case
%% candidate lambda's in simulation
e=(log(100)-log(.0001))/99;
lambda=exp(log(.0001):e:log(100)); 
%% number of folds in cross validation
nfold=10;
%% value of alpha
alpha1=1;   % For Lasso, alpha=1.
%% optimal lambda and alpha obtained vrom cross validation
[lambda_opt alpha_opt r]=cv_logitlap( Xm,y,L,alpha1,lambda,nfold);  
%% theta_hat according to logistic regression using optimal
%% values of lambda and alpha
 theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save logisdataLasso1 % save all work space data 
 bar(theta_hat); % bar plot the theta_hat values
xlabel('The dimension N');
title('Logistic Regression');
[Xm1,y1,L1]=SimLassoData(beta1,200,100); % testing data
lerror=logisSimError(Xm1,y1,theta_hat); % test the errors

