%% This script is to simulate the first Lasso test 50 times and get the error 
clear all; 
clc
load PPIdata; % data file containing X, y, Adjacency matrix A
lambda_opt=.2900410848; %% optimal lambda got by CV
alpha1=.2;
logisError=zeros(10,1);
Sp=zeros(10,1);
Sen=zeros(10,1);
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
    miu1(i,:)=miu;    sig1(i,:)=sig;
end
X1=(X'-miu1)./sig1;  % To normalize the matrix X
X1=X1';  xnot1=X(:,9:83); xnot10=X(:,1:72);
 %% Set the rest xnot's
 for i=2:(10-1)
    fnw=['xnot',int2str(i),'=[X(:,1:(floor(83/10))*(',int2str(i),'-1))  X(:, (floor(83/10)*',int2str(i),'+1):83)];'];
    eval(fnw);
 end
  %% Set ynot in a similar way
 ynot1=y(9:83,1);
   ynot10=y(1:72,1);

 for i=2:(10-1)
    fnw=['ynot',int2str(i),'=[y(1:(floor(83/10))*(',int2str(i),'-1)); y((floor(83/10)*i+1):83, 1)];'];
    eval(fnw);
 end
 for i=1:10
    fnw=['theta_hat=LogitisLap(xnot',int2str(i),',ynot',int2str(i),',L,lambda_opt*alpha1,lambda_opt*(1-alpha1));'];
    eval(fnw);
    y_hat=(1./(1+exp(-transpose(theta_hat)*X))>.5);
    [Sp(i) Sen(i)]=GetSpSenY(y',y_hat); % sensitivity and specificity computation
    logisError(i)=mse(transpose(y)-y_hat);
end
%% PMSE and standard errors 
squaremean=mean(logisError);
a1=mean(logisError);
fprintf('The PMSE is: %15.6f \n',squaremean)
standerr=sqrt(sum((logisError-a1).^2))/sqrt(10);
fprintf('The standard error of PMSE is: %15.6f \n',standerr)
%% sensitivity and standard error 
sensitivity=mean(Sen);
fprintf('The sensitivity is: %15.6f \n',sensitivity)
standerrsen=sqrt(sum((Sen-sensitivity).^2))/sqrt(10);
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)
%% specificity and standard error
specificity=mean(Sp);
fprintf('The specificity is: %15.6f \n',specificity)
standerrsp=sqrt(sum((Sp-specificity).^2))/sqrt(10);
fprintf('The standard error is of spacificity: %15.6f \n',standerrsp)
save PMSEReal % save work space data
