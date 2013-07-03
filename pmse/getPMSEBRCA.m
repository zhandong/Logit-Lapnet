%% script to apply the algorithm to BRCA model and get squared mean error 
clear all; 
clc
load dataAppBrca; % data file that contains related parameters
alpha1=.2; % optimal alpha is 0.2
logisError=zeros(10,1); % array to store errors
Sp=zeros(10,1); % array to store specificity
Sen=zeros(10,1); % array to store sensitivity
D=diag(sum(A,2)); 
D1=sqrt(D);
[dimd,~]=size(D1);
for i=1:dimd
    if ( D1(i,i)~= 0 )
        D1(i,i)=1/D1(i,i);
    end
end
%% graphical Laplacian matrix correspondingly
L=D-A;  % Laplacian graphical matrix
L1=D1*L*D1;  % normalized laplacian matrix
[n p]=size(X'); 
miu=mean(X'); %% get the mean 
sig=var(X'); % variance 
sig=sqrt(sig); % square root of variance
sig1=zeros(n,p); 
miu1=zeros(n,p);
for i=1:n
    miu1(i,:)=miu;
    sig1(i,:)=sig;
end
X1=(X'-miu1)./sig1;  % normalize matrix X
X=X1'; 
xnot1=X(:,9:83);
xnot10=X(:,1:72);
 %% the rest xnot's
 for i=2:(10-1)
    fnw=['xnot',int2str(i),'=[X(:,1:(floor(83/10))*(',int2str(i),'-1))  X(:, (floor(83/10)*',int2str(i),'+1):83)];'];
    eval(fnw);
 end
  %% ynot obtaining 
 ynot1=y(9:83,1); 
 ynot10=y(1:72,1);
 for i=2:(10-1)
    fnw=['ynot',int2str(i),'=[y(1:(floor(83/10))*(',int2str(i),'-1)); y((floor(83/10)*i+1):83, 1)];'];
    eval(fnw);
 end
 for i=1:10
    fnw=['theta_hat=LogitisLap(xnot',int2str(i),',ynot',int2str(i),',L,lambda_opt*alpha1,lambda_opt*(1-alpha1));'];
    eval(fnw);
    y_hat=(1./(1+exp(-transpose(theta_hat)*X))>.5); % estimated responses
    [Sp(i) Sen(i)]=GetSpSenY(y',y_hat);
    logisError(i)=mse(transpose(y)-y_hat);
end
%% PMSE and standard errors 
squaremean=mean(logisError);
a1=mean(logisError);
fprintf('The PMSE is: %15.6f \n',squaremean)
standerr=sqrt(sum((logisError-a1).^2))/sqrt(10);
fprintf('The standard error of PMSE is: %15.6f \n',standerr)
%% sensitivity and standard penalized mean error
sensitivity=mean(Sen);
fprintf('The sensitivity is: %15.6f \n',sensitivity)
standerrsen=sqrt(sum((Sen-sensitivity).^2))/sqrt(10);% standard square root of sensitivity
fprintf('The standard error of sensitivity is: %15.6f \n',standerrsen)
%% specificity and standard error 
specificity=mean(Sp);
fprintf('The specificity is: %15.6f \n',specificity)
standerrsp=sqrt(sum((Sp-specificity).^2))/sqrt(10); % mean square root of specificity
fprintf('The standard error is of spacificity: %15.6f \n',standerrsp)
save PMSEAppBRCA 
