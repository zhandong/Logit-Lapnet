function [Xm,y,L]=SimLassoData(theta,n,p)
%% Generate necessary data in simulation of Lasso 
%% The inputs are theta, and the dimension of transcription factor: n by p
%% The outputs are the matrix Xm, the vector y, the Laplacian matrix L
Xm=zeros(n+10*n,p);  % Caiyan Li's paper: n=200, p=100
x=randn(n,p);
%% x is the matrix containing transcription factors which follow normal distributions
%% Xm is the matrix that contains all the genes that include the transcription
%% factors and the genes that the TFs regulate
%% i.e., x[i,j]--> [x[i,j]; randn(10,1)*sqrt(0.5)+0.7*x[i,j] ]
%% ++++++++++++++++++++++++++++++++++++++++++ get expression matrix Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j);
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];
    end
end
%% ++++++++++++++++++++++++++++++++++++++++++ get y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % logistic regression
%% ++++++++++++++++++++++++++++++++++++++++++
%% the graphical adjacency matrix 
%% ++++++++++++++++++++++++++++++++++++++++++ get A
 A=zeros(n+10*n,n+10*n);
 for i=1:11:(n+10*n-10)
     A((i+1):(i+10),i)=ones(10,1);   
     A(i,(i+1):(i+10))=ones(1,10);
 end
 %%++++++++++++++++++++++++++++++++++++++++++++++++++++
%% degree matrix of A, i.e., D
D=diag(sum(A,2));
%% graphical Laplacian matrix 
L=D-A; 
