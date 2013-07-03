clear all; 
clc;
load RealData; %% load data generated from wet lab
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
L1=D1*L*D1;  % normalize laplacian matrix
[n p]=size(X'); 
miu=mean(X'); % get mean of the real data
sig=var(X'); % variance of the data
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p);
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;   
    sig1(i,:)=sig;
end
X1=log(1+X'); % log the data X
X1=(X1-miu1)./sig1; % normalize the matrix X for use
alpha1=0.2;
%% get least upper bound of lambda
lammax=getLambMax(X1,y,alpha1); 
e=(log(lammax)-log(1))/19; 
lambda=exp(log(1):e:log(lammax)); 
%% +++++++++++++++++ the number of folds in cross validation
nfold=10;
%% value of alpha's
%% get optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% get theta_hat according to logistic regression using optimal
%% values for lambda and alpha
 theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save AppLogNorm  % save work space data
 bar(theta_hat); % bar plot estimated theta
xlabel('The dimension N');
title('Logistic Regression')
lerror=logisSimError(X,y,theta_hat); % get errors of coefficients


