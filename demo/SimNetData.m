function [Xm,y,L]=SimNetData(theta,n,p)
%% Generate data in simulations for Logit-Lap Net 
%% Inputs are theta, dimension of transcription factors
%% Outputs are gene expression matrix Xm, response vector y, Laplacian matrix L
Xm=zeros(n+10*n,p);  %% Hongzhe Li's paper: n=200, p=100
x=randn(n,p);
%% Matrix x contains transcription factors that follow normal distribution
%% Matrix Xm contains the genes that include transcription factors and the 
%% genes that each of the TF regulates
%% i.e., x[i,j]-->[x[i,j]; randn(10,1)*sqrt(0.5)+0.7*x[i,j] ]
%% First, get Xm and y
%% +++++++++++++++++++++++++++++++++++++++++++++++  Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j); % each TF regulates 10 genes
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];
    end
end
%% +++++++++++++++++++++++++++++++++++++++++++++++  y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % for logistic regression model
%% +++++++++++++++++++++++++++++++++++++++++++++++
%% the graphical adjacency matrix 
%% +++++++++++++++++++++++++++++++++++++++++++++++ 
 A=zeros(n+10*n,n+10*n); % each TF regulates 10 genes
 for i=1:11:(n+10*n-10)
     A((i+1):(i+10),i)=ones(10,1);
     A(i,(i+1):(i+10))=ones(1,10);
 end
%% +++++++++++++++++++++++++++++++++++++++++++++++
%% degree matrix of A:  D
D=diag(sum(A,2));
%% graphical Laplacian matrix 
L=D-A; 