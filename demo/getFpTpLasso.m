%% script to simulate the four experiments for Lasso and get ROC( Receiver
%% operator Curve)
clear all;
load coeff1; % coefficients in four simulations
% load coeff2; % coefficients in Simulation 1 
% load coeff3; % coefficients in Simulation 1 
% load coeff4; % coefficients in Simulation 1 
alpha1=1; % for Lasso, alpha is 1
[Xm,y,L]=SimLassoData(beta1,200,100); % get simulation related data
lammax=getLambMax(Xm',y,alpha1);  % get maximal lambda
e=(log(lammax)-log(.0001))/19;
% set lambda candidates
lambda=exp(log(.0001):e:log(lammax)); 
Fpr=zeros(50,1); % store false positive
Tpr=zeros(50,1); % store true positive
for i=1:50
     % training handle
     fnw2=['theta_hat=LogitisLap(Xm,y,L,lambda(',int2str(i),')*alpha1,lambda(',int2str(i),')*(1-alpha1));'];
     eval(fnw2);
     [Fpr(i) Tpr(i)]=GetFPTP(beta1,theta_hat, 0.0001); % FP and TP computing
end
save FPTPLasso % save work space data

