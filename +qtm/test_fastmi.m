load GSM3308547_GSM330854.mat
%y = sce.list_cell_attributes{2};
%X = sc_transform(sce.X);
[~, X, g] = sc_splinefit(sce.X, sce.g);
%X = sc_norm(X);
X = sc_transform(X,"type","PearsonResiduals");
% X = sc_impute(X,"type","MAGIC");
n = 2000;
X = X(1:n,:);
g = g(1:n);
y = X(18,:);

X(18,:)=[];
g(18)
g(18)=[];

% idx=randperm(length(g));
% g=g(idx);
% X=X(idx,:);

h = 1.2;

K = 15;

data = [X; y];
R0  = qtm.FastPairMI(data, h);

R = R0(1:end-1,1:end-1)/(K-1);
J = R0(end,1:end-1);

tic;
fun = @(alpha)howmany(alpha,R,J) - K;
alphasol = fzero(fun,[0 1]);
% alphasol = 0.3;
[~,xsol] = howmany(alphasol,R,J);
g(find(xsol.BestX))
toc;

tic;
B = lasso(X',y');
g(B(:,80)~=0)
toc;

intersect(g(find(xsol.BestX)), g(B(:,80)~=0))


function [n,result] = howmany(alpha,R,J)
    Q = qubo((1-alpha)*R - alpha*diag(J));
    result = solve(Q);
    n = sum(result.BestX);
end


