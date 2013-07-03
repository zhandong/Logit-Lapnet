clear all; 
clc;
%% data from wet lab
load RealData; 
A=Afil; 
X=Xfil; 
y=yfil;
D=diag(sum(A,2)); 
D1=sqrt(D); 
[dimd,~]=size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
%% graphical Laplacian matrix 
L=D-A;  
L1=D1*L*D1;  % normalized laplacian matrix
[n p]=size(X'); 
miu=mean(X'); % mean of real data
sig=var(X'); % variance of the data
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p); 
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;   
    sig1(i,:)=sig;
end
X1=log(1+X'); % log without norm
%% get the value of alpha's
alpha1=1; % Lasso case
%% get least upper bound of lambda
lammax=getLambMax(X',y,alpha1); 
e=(log(lammax)-log(1))/19;
lambda=exp(log(1):e:log(lammax)); 
%% ++++++++++++++++++ get the number of folds in cross validation
nfold=10;
%% get the optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% get the estimated theta using logistic regression via optimal
%% values for lambda and alpha
 theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save LogWoNormLasso; % save work space data 
 bar(theta_hat);
 xlabel('The dimension N'); 
title('Logistic Regression');
lerror=logisSimError(X,y,theta_hat);


