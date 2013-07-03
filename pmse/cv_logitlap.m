function [lambda_opt alpha_opt r]=cv_logitlap(X,y,L,alpha,lambda,nfold)
%% function using cross validation (CV) to get optimal lambda and alpha,
%% which is the so called generative process.  
%% nfold: the number of folds in CV
%% X: gene expression matrix, y: response vector; L: graphical Laplacian matrix
%% alpha and lambda: variables in loss function that to be determined, refer
%% to equation (3) for information
dimG=size(X,1); % get dimension G 
dimn=size(y,1); % get dimension n
%% get number of candidate lambda's
[m,n]=size(lambda);
if m>n
    sizelam=m;
else
    sizelam=n;
end
%% get length of candidate alpha's. 
[m,n]=size(alpha);
if m>n
    sizealpha=m;
else
    sizealpha=n;
end
%% divide X and y into several folds
 for i=1:nfold
      fnw=['x',int2str(i),'=X(:,(dimn/nfold*',int2str(i),'-dimn/nfold+1):(dimn/nfold*',int2str(i),'));'];
      eval(fnw);
      fnw=['y',int2str(i),'=y((dimn/nfold*',int2str(i),'-dimn/nfold+1):(dimn/nfold*',int2str(i),'),1);'];
      eval(fnw);
 end
 %% get head and rear xnot first
 xnot1=X(:,dimn/nfold+1:dimn); 
 fnw=['xnot',int2str(nfold),'=[X(:,1:dimn-dimn/nfold)];'];
 eval(fnw);
 %% get xnot's in the middle of the matrix
 for i=2:(nfold-1)
     fnw=['xnot',int2str(i),'=[X(:,1:(dimn/nfold)*(',int2str(i),'-1))  X(:, (dimn/nfold*',int2str(i),'+1):dimn)];'];
     eval(fnw);
 end
 %% get ynot also
 ynot1=y(dimn/nfold+1:dimn,1);
 fnw=['ynot',int2str(nfold),'=[y(1:dimn-dimn/nfold,1)];'];
 eval(fnw);
 for i=2:(nfold-1)
     fnw=['ynot',int2str(i),'=[y(1:(dimn/nfold)*(',int2str(i),'-1)); y((dimn/nfold*i+1):dimn, 1)];'];
     eval(fnw);
 end
%% the followings are cross validation process to obtain optimal alpha and lambda
for k=1:sizealpha
 alphacan=alpha(k);    
 for i=1:nfold
     for j=1:sizelam
         %% Theta=LogitisLap (xnot_i,ynot_i,L,lambda(j)*alphacan,lambda(j)*(1-alphacan));
         fnw=['Theta=LogitisLap (xnot',int2str(i),',ynot',int2str(i),', L,lambda(',int2str(j),')*alphacan,lambda(',int2str(j),')*(1-alphacan));'];
         eval(fnw);
         %% r: matrix that contains related residual=|| y_i-(1/exp(-theta'*X_i)>0.5) ||_2
         fnw=['r(',int2str(i),',',int2str(j),')=norm((transpose(y',int2str(i),')-(1./(1+exp(-transpose(Theta)*x',int2str(i),'))>.5)),2);'];
         eval(fnw);
     end
 end
 %% get minimal residual and the index among all candidates
 rmin(k)= max(mean(r,1)); % set maximal mean to rmax
 optindex=1;  % in case, the optimal lambda is the first one
 for j=1:sizelam
     if rmin>mean(r(:,j))
         rmin(k)=mean(r(:,j)); 
         optindex=j;  % get index of the minimal residual
     end
 end
 %% get optimal lambda
 lambdaopt(k)=lambda(optindex);  % get optimal lambda 
end
[err Ind]=min(rmin); 
alpha_opt=alpha(Ind); 
lambda_opt=lambdaopt(Ind);

    