%% codes for application of Logit-Lapnet algorithm. Version 1.1.1, developed
%% during Mar. to May. 2013 by Zhandong Liu's lab. For detailed information,
%% refer to the paper by Zhang, W., Wan, Y. W., Allen, G., Liu, Z.D. et al, 
%% 'Molecular pathway identification using biological network-regularized 
%% logistic models'. International Conference on Intelligent Biology and Medicine
%% (ICIBM 2013). Also to appear in a special issue of BMC Genomics

clear all; 
clc;
load RealData; % data from the file that contains X, y, Adjacency A
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
miu=mean(X'); % mean of real data
sig=var(X'); % variance of real data
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p); 
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;   
    sig1(i,:)=sig;
end
X1=log(1+X'); % log without norm
alpha1=0.2;
%% get the least upper bound of lambda, here to X is smaller
lammax=getLambMax(X',y,alpha1); 
e=(log(lammax)-log(1))/19;
lambda=exp(log(1):e:log(lammax)); 
%% ++++++++++++++++++ Set the number of folds in cross validation
nfold=10;
%% give the value of alpha's
%% get optimal lambda and alpha through cross validation 
[lambda_opt alpha_opt r]=cv_RealApplogitlap( X1',y,L,alpha1,lambda,nfold);
%% get estimated coefficients according to logistic regression using optimal
%% values for lambda and alpha
 theta_hat=LogitisLap(X1',y,L,lambda_opt*alpha_opt,lambda_opt*(1-alpha_opt));
 save AppLogWoNorm % save work space data
 bar(theta_hat); % bar plot the estimated coefficient
xlabel('The dimension N'); 
title('Logistic Regression')
lerror=logisSimError(X,y,theta_hat);


