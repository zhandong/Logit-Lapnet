function theta = LogitisLap (X,y,L,lambda1,lambda2)
%% get the argmin of loss function -sum_i{y_i*log(g(theta'*x_i))+(1-y_i)*log(1-g(theta'*x_i))}
%% +lambda1 |theta|_1+lambda2 theta' A theta
%% Inputs are: y - an n by 1 vector; X - a G by n matrix
%% L is the graph Laplacian matrix
%% lambda1 and lambda2 are parameters known in the loss function  
%% The output theta is the optimal solution of the covex optimization problem
%% The dimension of theta: G by 1
[dimn,~]=size(y); % the number of dimension dimn
[dimG,m]=size(X); % the number of dimension dimG, m is also useful in the following
%%  To use cvx to solve the convex problem
cvx_begin
    variables theta(dimG) ;
    minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( [zeros(1,m); -theta'*X]) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
cvx_end
    
    
