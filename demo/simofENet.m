clear all; 
clc
load coeff1; % coefficients in Simulation 1
% load coeff2; % coefficients in Simulation 2 
% load coeff3; % coefficients in Simulation 3 
% load coeff4; % coefficients in Simulation 4 
%% get gene expression matrix Xm, response vector y, graphical Laplacian matrix L 
[Xm,y,L]=SimEnetData(beta1,200,100); % LogitLap net simulate data
%% get candidate lambda's in simulation
e=(log(100)-log(.0001))/99;
lambda=exp(log(.0001):e:log(100)); 
%% ++++++++++++++++++ number of folds in cross validation
nfold=10;
%% optimal values of alpha 
alpha1=0.2;   % alpha here can also be a vector that contains candidates  
%% obtain optimal lambda and alpha via CV, r is the resudial
%% inputs are gene expression matrix Xm, response vector y; graphical Laplacian
%% L, candidates alpha, lambda and the # of nfold
[lambda_opt alpha_opt r]=cv_logitlap(Xm,y,L,alpha1,lambda,nfold);  
%% estimated coefficients from logistic regression 
%% please see equation (3) of the paper for parameters in loss function
theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdataNet2 % save work space data for reusage
bar(theta_hat); % bar plot estimated coefficients
xlabel('The dimension N');
title('Logistic Regression')
[Xm1,y1,L1]=SimNetData(beta1,200,100);
lerror=logisSimError(Xm1,y1,theta_hat); % the L2 norm of the errors

