clear all;
load coeff1; % different coefficients in 4 simulations
% load coeff2;  % the 2nd simulation
%  load coeff3; % the 3rd simulation
%  load coeff4; % the 4th simulation
%% Xm, y, A in simulation
[Xm,y,L]=SimLassoData(beta1,200,100); % lasso case
%% candidate lambda's in simulation
 e=(log(100)-log(.0001))/99;
 lambda=exp(log(.0001):e:log(100)); 
%% ++++++++++++++++++ the number of folds
nfold=10;
%% get value of alpha, for lasso is 1 
alpha1=1; 
%% get optimal lambda through cross validation process
[lambda_opt r]=cv_logitlap( Xm,y,L,alpha1,lambda,nfold);  
%% get estimated theta_hat via logistic regression
 theta_hat=LogitisLap(Xm,y,L,lambda_opt*alpha1,lambda_opt*(1-alpha1));
 save logisdata 
%% Plot first 100 elements 
plot([1:100] , theta_hat(1:100),'o',[1:100],beta1(1:100),'*')
xlabel('The dimension N');
title('Logistic Regression');
legend('Thetahat','Beta','Location','NorthEast');
[Xm1,y1,L1]=LogisSimData(beta1,200,100);
lerror=logisSimError(Xm1,y1,theta_hat);



