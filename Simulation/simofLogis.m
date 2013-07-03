%% The function of simulations
%% The updated files are simofLogisNet, simofLogisENet, simofLogisLasso 
clear all; 
clc;
load coeff1; % different coefficients in 4 simulations
% load coeff2;  % the 2nd simulation
%  load coeff3; % the 3rd simulation
%  load coeff4; % the 4th simulation

%% to get Xm, y, A in simulation
[Xm,y,L]=SimEnetData(beta1,200,100); % Elastic net
%% set candidate lambda's in simulation
e=(log(100)-log(.0001))/2;
lambda=exp(log(.0001):e:log(100)); 
%% ++++++++++++++++++ Give number of folds that divided to, theoretically any integer is okay
nfold=2;
%% give value of candidate (optimal) alpha 
alpha1=1;  %% for Lasso, alpha is 1
%% to get the optimal lambda through cross validation
[lambda_opt alpha_opt r]=cv_logitlap( Xm,y,L,alpha1,lambda,nfold);  
% get the theta_hat according to logistic regression
theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdata 
%% Plot first 100 elements 
plot([1:100] , theta_hat(1:100),'o',[1:100],beta1(1:100),'*')
xlabel('The dimension N');
title('Logistic Regression');
legend('Thetahat','Beta','Location','NorthEast');
[Xm1,y1,L1]=SimEnetData(beta1,200,100);
lerror=logisSimError(Xm1,y1,theta_hat);

