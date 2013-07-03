clear all; 
clc
load coeff1; % load coefficients in simulation1
% load coeff2;   % coefficients in simulation 2
% load coeff3;   % coefficients in simulation 3
% load coeff4;   % coefficients in simulation 4
%% get Xm, y, Laplacian matrix in simulation
[Xm,y,L]=SimNetData(beta1,200,100); % simulate Logit-lap net
%% candidate lambda's in simulations
e=(log(100)-log(.0001))/99;
lambda=exp(log(.0001):e:log(100)); 
%% number of folds in cross validation, ten-fold CV for instance
nfold=10;
%% set candidate values of alpha 
alpha1=0.2:.2:.6;  
%% get optimal lambda and alpha through cross validation
[lambda_opt alpha_opt r]=cv_logitlap( Xm,y,L,alpha1,lambda,nfold);  
%% estimated coefficient according to logistic regression by optimal alpha and lambda
theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdataNet 
%% barplot theta_hat
bar(theta_hat);
xlabel('The dimension N');
title('Logistic Regression');
[Xm1,y1,L1]=SimNetData(beta1,200,100);
lerror=logisSimError(Xm1,y1,theta_hat);

