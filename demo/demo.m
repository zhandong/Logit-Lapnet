%% demo file is to test the installation of your CVX package as well as
%% proper usage of our package. In command interface of Matlab, type demo
%% to see what will display. the demo file runs lasso actually.
clear all; 
clc
load coeff1; % simulation1 for example
%% get gene expression matrix Xm, response vector y, graphical Laplacian matrix L 
[Xm,y,L]=SimLassoData(beta1,200,100); % lasso case
%% candidate lambda's in demo example
e=(log(100)-log(.0001))/99;
lambda=exp(log(.0001):e:log(100)); 
%% ++++++++++++++++++ ten-fold cross validation
nfold=10;
%% optimal values of alpha 
alpha1=.3;   % alpha here can be a vector also
%% obtain optimal lambda and alpha 
%% inputs are gene expression matrix Xm, response vector y; graphical Laplacian
%% L, candidates alpha, lambda and the number of folds that divided
disp('--- Demo program, wait to enjoy the wordy interface ....');
fprintf('\n')
disp('--- Wordy? Sure, you installed CVX successfully ....');
fprintf('\n')
[lambda_opt alpha_opt r]=cv_logitlap(Xm,y,L,alpha1,lambda,nfold);  
%% estimated coefficients from logistic regression 
theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdataDemo % reusable data
bar(theta_hat); % bar plot estimated coefficients
xlabel('The dimension N');
title('Logistic Regression')
[Xm1,y1,L1]=SimNetData(beta1,200,100);
lerror=logisSimError(Xm1,y1,theta_hat); % L2 norm of the errors

