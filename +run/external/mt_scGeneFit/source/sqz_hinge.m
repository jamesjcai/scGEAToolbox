%Filename:      sqz_hinge.m
%Last edit:     Feb 11 2019
%Description:   Implements the LP to select markers that preserve 
%               classification structure using Matlab's linprog'
%Inputs:
%               -Delta: 
%     
%               a d x N array of constraints, d is the dimension
%               of the space and N is the number of constraints
% 
%               -epsilon: 
% 
%               a hinge parameter
%
%               -dim_target: 
%     
%               number of markers to select
%
%               -lambda: 
%     
%               a N x d array stating a penalization for each constraint 
%               (typically all ones or inversely proportional to the number of 
%               on the number of points in the cluster) 
%Outputs:
%               -M: 
% 
%               a dim x 1 binary array. M(t)=1 implies that t is a selected marker
%
%               
% 
%Documentation:
% 
% -------------------------------------------------------------------------
function M=sqz_hinge(Delta, epsilon, dim_target,lambda)
[d,N]=size(Delta);

f=horzcat(zeros(1,d), ones(1,N));
lb=zeros(1,d+N);

b=-epsilon*ones(N,1)./lambda;
b=vertcat(b, dim_target);

aux=-Delta';
A=horzcat(aux, -speye(N,N));
A=vertcat(A, [ones(1,d), zeros(1,N)]);

ub=[ones(1,d), inf*ones(1,N)];
opts=optimoptions('linprog', 'Display', 'iter', 'OptimalityTolerance', 1e-4,'ConstraintTolerance', 1e-3);

M=linprog(f,A,b,[],[],lb,ub, opts);
M=M(1:d);
end
