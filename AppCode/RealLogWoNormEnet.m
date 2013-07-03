clear all; 
clc;
load RealData; % data file from wet lab
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
[m,~]=size(A); 
L=eye(m); % Elastic: L is identity
[n p]=size(X'); 
miu=mean(X'); % get mean of data
sig=var(X'); % variance of data
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p); 
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;   
    sig1(i,:)=sig;
end
X1=log(1+X'); % log without norm
alpha1=0.2;
%% get least upper bound of lambda, here use X whose is smaller
lammax=getLambMax(X',y,alpha1); 
e=(log(lammax)-log(1))/19;
lambda=exp(log(1):e:log(lammax)); 
%% ++++++++++++++++++ Set the number of folds in cross validation
nfold=10;
%% get value of alpha's
%% get the optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% get the estimated theta via the logistic regression using optimal
%% values of lambda and alpha
 theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save LogWoNormEnet;  % save all the variables for reusage
 bar(theta_hat);
 xlabel('The dimension N');
title('Logistic Regression');
lerror=logisSimError(X,y,theta_hat);


