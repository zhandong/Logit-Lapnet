%% script to simulate experimentd for Elastic net so as to get the ROC
clear all;
load coeff1; % coefficients in four simulations
% load coeff2; % coefficients in Simulation 2 
% load coeff3; % coefficients in Simulation 3 
% load coeff4; % coefficients in Simulation 4 

 alpha1=.2; % via CV, choose alpha1 as .2
 [Xm,y,L]=SimEnetData(beta1,200,100);
 lammax=getLambMax(Xm',y,alpha1); % lambda_max from equation (6) in the paper
 e=(log(lammax)-log(1))/19;
 %% candidate lambda set
lambda=exp(log(1):e:log(lammax)); 
Fpr=zeros(20,1);
Tpr=zeros(20,1);
for i=1:20
     fnw2=['theta_hat=LogitisLap(Xm,y,L,lambda(',int2str(i),')*alpha1,lambda(',int2str(i),')*(1-alpha1));'];
     eval(fnw2);
    [Fpr(i) Tpr(i)]=GetFPTP(beta1,theta_hat, 0.0001); % sensitivity and specificity computation
end
save FPTPENet %% Fpr and Tpr are stored so that ROC could be obtained then

