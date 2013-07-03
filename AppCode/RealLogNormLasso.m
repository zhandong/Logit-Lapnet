clear all; 
clc;
load RealData; 
%% data file that contains X, y, adjacency matrix A
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
%% The graph Laplacian matrix correspondingly
L=D-A;  
L1=D1*L*D1;  % get normalized laplacian matrix
[n p]=size(X'); 
miu=mean(X'); % get the mean of real data
sig=var(X'); % variance of the real data
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p); 
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;   
    sig1(i,:)=sig;
end
X1=log(1+X'); % log without norm
X1=(X1-miu1)./sig1; % normalize the matrix X
alpha1=1; % Lasso case
%% get the least upper bound of lambda, to X is smaller here
lammax=getLambMax(X',y,alpha1); 
e=(log(lammax)-log(1))/19;
lambda=exp(log(1):e:log(lammax)); 
%% ++++++++++++++++++ Set the number of folds in cross validation
nfold=10;
%% provide the value of alpha's
%% get the optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% To get the theta_hat according to the logistic regression using optimal
%% values for lambda and alpha
 theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save LogNormLasso;  
 bar(theta_hat);
 xlabel('The dimension N');
title('Logistic Regression');
lerror=logisSimError(X,y,theta_hat);


