function [Xm,y,L]=SimNetData(theta,n,p)
%% Generate necessary data in the simulation for Logit-Lap Net 
%% Inputs are theta, the dimension of transcription factors
%% The outputs are the matrix Xm, the vector y, the Laplacian matrix L
Xm=zeros(n+10*n,p);  %% In Li's paper, n=200, p=100
x=randn(n,p);
%% Matrix x contains transcription factors which follow normal distribution
%% Matrix Xm contains all the genes that include the transcription
%% factors and the genes that TFs regulate

%% get the Xm and y
%% +++++++++++++++++++++++++++++++++++++++++++++++++ Set Xm
for i=1:n
    for j=1:p
        tmp=randn(10,1)*sqrt(.5)+.7*x(i,j); % each TF regulates 10 genes
        Xm((11*i-10):(11*i),j)=[x(i,j);tmp];% following Caiyan Li's paper
    end
end
%% +++++++++++++++++++++++++++++++++++++++++++++++++ get y 
y=(1./(1+exp(-theta'*Xm)))'>rand(1,p)'; % logistic model
%% +++++++++++++++++++++++++++++++++++++++++++++++++
%% Get the graph adjacency matrix 
%% +++++++++++++++++++++++++++++++++++++++++++++++++ adjacency matrix A
 A=zeros(n+10*n,n+10*n); % each TF regulates 10 genes
 for i=1:11:(n+10*n-10)
     A((i+1):(i+10),i)=ones(10,1);
     A(i,(i+1):(i+10))=ones(1,10);
 end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++
%% get degree matrix of A and set it to D
D=diag(sum(A,2));
%% graphical Laplacian matrix 
L=D-A; 