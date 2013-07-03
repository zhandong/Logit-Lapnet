%% The script for real application running: Logistic model for real network regularization
clear all; 
clc
load PPIdata; % data contains X, y, and adjacency matrix A
y=y(:,3)-1; % response vector
D=diag(sum(A,2));  % degree matrix of simulation
D1=sqrt(D);  % square root of degree-normalize
[dimd,~]=size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
%% graphical Laplacian matrix
L=D-A; 
L1=D1*L*D1;  % normalize laplacian matrix
[n p]=size(X');
miu=mean(X');   % mean of real data
sig=var(X');    % variance of the real data
sig=sqrt(sig);  % square root of variance
sig1=zeros(n,p);
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;
    sig1(i,:)=sig;
end
X1=(X'-miu1)./sig1;  %  normalize the gene expression matrix X
%% set candidate lambda's in application
alpha1=0.2;
lammax=getLambMax(X1,y,alpha1);  % upper bound of lambda
e=(log(lammax)-log(1))/19;  % set lambda_min 
lambda=exp(log(1):e:log(lammax)); 
%% ++++++++++++++++++ the number of folds in cross validation
nfold=10;
%% give the value of alpha's
%% get the optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% get the theta_hat according to logistic regression using optimal
%% lambda and alpha's
theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
save logisdataApp % save workspace data
bar(theta_hat); 
xlabel('The dimension N');
title('Logistic Regression'); 
lerror=logisSimError(X,y,theta_hat);

