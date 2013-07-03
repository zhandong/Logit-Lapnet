function [Xm,y,L]=SimEnetData(theta,n,p)
%% Generate testing and training data in simulations of Elastic Net 
%% The inputs are theta, the dimensions of TF ( n by p)
%% The outputs are expression Xm, response y, graphical Laplacian matrix L
Xm=zeros(n+10*n,p);  
x=randn(n,p);
%% x is the matrix that contains transcription factors which follow normal distributions
%% Xm is the matrix that contains all the genes including the transcription
%% factors and the genes that TFs regulate
%% i.e., x[i,j]--> [x[i,j]; randn(10,1)*sqrt(0.5)+0.7*x[i,j] ]
%% First, get the Xm and y
%% +++++++++++++++++++++++++++++++++++++++++++++++++ get Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j);
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];
    end
end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++ get y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % logistic model
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++
%% The graphical Laplacian matrix of Elastic net is identity
L=eye(n+10*n,n+10*n);