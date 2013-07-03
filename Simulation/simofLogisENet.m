clear all;
clc
load coeff1; % the 1st simulation
% load coeff2;  % the 2nd simulation
% load coeff3; % the 3rd simulation
% load coeff4; % the 4th simulation
%% Xm, y, L in simulations
[Xm,y,L]=SimEnetData(beta1,200,100); % simulate Elastic net
%% candidate lambda's in simulations
e=(log(100)-log(.0001))/99;
lambda=exp(log(.0001):e:log(100)); 
%% ++++++++++++++++++ number of folds that data is divided in cross validation
nfold=10;
%% get value of candidate alpha's 
alpha1=0.2:.2:.6; 
%% get the optimal lambda and alpha through cross validation
[lambda_opt alpha_opt r]=cv_logitlap( Xm,y,L,alpha1,lambda,nfold);  
%% get the theta_hat according to the logistic regression using optimal
%% values for lambda and alpha
 theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save logisdataENet; % work space data
 bar(theta_hat); % bar plot the estimated coefficient
xlabel('The dimension N');
title('Logistic Regression');
[Xm1,y1,L1]=SimEnetData(beta1,200,100); %% testing data
lerror=logisSimError(Xm1,y1,theta_hat);

