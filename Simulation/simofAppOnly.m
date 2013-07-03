%% This is the script for real application running---Logistic model for real network regularization
clear all; 
clc
load PPIOnlydata; % data generated from the file which contains X, y, Adjacency A
D=diag(sum(A,2));
D1=sqrt(D);
[dimd,~]=size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
%% The graph Laplacian matrix correspondingly
L=D-A; 
L1=D1*L*D1;  % normalize the Laplacian
[n p]=size(X');
miu=mean(X'); % get mean of real data
sig=var(X'); % variance of real application data
sig=sqrt(sig); % square root of variances
sig1=zeros(n,p);
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;
    sig1(i,:)=sig;
end
X1=(X'-miu1)./sig1;  % normalize the matrix X
%% set candidate lambda's in application
alpha1=0.2;
lammax=getLambMax(X1,y,alpha1); %% to get the upper bound of lambda
e=(log(lammax)-log(.001))/19;
lambda=exp(log(.001):e:log(lammax)); 
%% ++++++++++++++++++ Set number of folds in CV
nfold=10;
%% to get the optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L1,alpha1,lambda,nfold);
%% get estimated theta_hat according to logistic regression using optimal
%% values for lambda and alpha, which finish the estimation process
theta_hat=LogitisLap(X1',y,L1,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdataApp % work space data for reusage
bar(theta_hat);
xlabel('The dimension N');
title('Logistic Regression')
lerror=logisSimError(X,y,theta_hat); % standard errors or coefficients


