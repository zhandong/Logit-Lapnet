%% script to simulate the four experiments for Elastic and get ROC ( Receiver
%% operator curve)
clear all;
load coeff1; % coefficients in four simulations
% load coeff2; % coefficients in Simulation 1 
% load coeff3; % coefficients in Simulation 1 
% load coeff4; % coefficients in Simulation 1 
alpha1=.2; % via CV, optimal alpha is .2
[Xm,y,L]=SimNetData(beta1,200,100);
lammax=getLambMax(Xm',y,alpha1); 
e=(log(lammax)-log(1))/19;
lambda=exp(log(1):e:log(lammax)); 
Fpr=zeros(50,1); % false positive
Tpr=zeros(50,1); % true positive
for i=1:50
    %% for different testing data, use fnw as handle
     fnw2=['theta_hat=LogitisLap(Xm,y,L,lambda(',int2str(i),')*alpha1,lambda(',int2str(i),')*(1-alpha1));'];
     eval(fnw2); % evaluate the handle
    [Fpr(i) Tpr(i)]=GetFPTP(beta1,theta_hat, 0.0001); % FP and TP computing
end
save FPTPLapNet % work space data for reusage

