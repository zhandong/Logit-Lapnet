function [Xm,y,L]=SimEnetData(theta,n,p)
%% generate necessary data in simulations for Elastic Net 
%% The inputs are theta, the dimensions of TF (n by p)
%% The outputs are matrix Xm, vector y, graphical Laplacian matrix L
Xm=zeros(n+10*n,p);  
x=randn(n,p);
%% x is the matrix that contains transcription factors following normal distribution
%% Xm is the gene expression matrix
%% i.e., x[i,j]--> [x[i,j]; randn(10,1)*sqrt(0.5)+0.7*x[i,j] ]
%% get Xm and y firstly
%% +++++++++++++++++++++++++++++++++++++++++++++++++ get Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j);
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];
    end
end
%% +++++++++++++++++++++++++++++++++++++++++++++++++ get y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % logistic regression
%% +++++++++++++++++++++++++++++++++++++++++++++++++
%% The graphical Laplacian matrix of Elastic net is identity
L=eye(n+10*n,n+10*n);