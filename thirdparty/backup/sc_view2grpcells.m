function [s,c]=sc_view2grpcells(X0,X1,s)

if size(X0,1)~=size(X1,1)
    error('X0 and X1 should contain same number of genes.');
end

if nargin<3
    s=sc_tsne([X0, X1],3);
end
n0=size(X0,2);
n1=size(X1,2);
c=[ones(n0,1);2*ones(n1,1)];
% i_gscatter3(s,c);
% i_myscatter(s,c);
scatter3(s(1:n0,1),s(1:n0,2),s(1:n0,3),10);
hold on
scatter3(s(n0+1:end,1),s(n0+1:end,2),s(n0+1:end,3),10);
% xlabel('t-SNE 1'); ylabel('t-SNE 2'); zlabel('t-SNE 3');

