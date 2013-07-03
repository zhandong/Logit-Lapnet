function [Xm,y,L]=SimLassoData(theta,n,p)
%% Generate testing and training data in simulation of Lasso case
%% The inputs are theta, the dimension of transcription factor, say, n by p
%% The outputs are the matrix Xm, the vector y, the Laplacian matrix L
Xm=zeros(n+10*n,p);  % following Hongzhe Li's work: set n=200, p=100
x=randn(n,p);
%% x is the matrix containing transcription factors 
%% Xm is the matrix that contains the genes including the transcription
%% factors and the genes that each TF regulates
%% i.e., x[i,j] --> [x[i,j]; randn(10,1)*sqrt(0.5)+0.7*x[i,j] ]
%% First, get Xm and y
%% ++++++++++++++++++++++++++++++++++++++++++ get expression matrix Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j);
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];
    end
end
%% ++++++++++++++++++++++++++++++++++++++++++  get response vector y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % for logistic regression
%% ++++++++++++++++++++++++++++++++++++++++++
%% Get graphical adjacency matrix 
%% ++++++++++++++++++++++++++++++++++++++++++ get adjacency matrix A
 A=zeros(n+10*n,n+10*n);
 for i=1:11:(n+10*n-10)
     A((i+1):(i+10),i)=ones(10,1);   
     A(i,(i+1):(i+10))=ones(1,10);
 end
 %%++++++++++++++++++++++++++++++++++++++++++++++++++++
%% degree matrix of A, which is set to D
D=diag(sum(A,2));
%% graphical Laplacian matrix 
L=D-A; 
