function [G]=sc_knngraph(s,k,plotit,method)
%Generate KNN group network from cell embeddings

if nargin<4, method=1; end
if nargin<3, plotit=false; end
if nargin<2 || isempty(k), k=4; end

switch method
    case 1
        [mIdx] = knnsearch(s,s,'K',k+1);
        Graph=mIdx';
    case 2
        pw1=fileparts(mfilename('fullpath'));
        pth=fullfile(pw1,'thirdparty/k-NN-code');
        addpath(pth);
        kneighbors=k; % number of neighbors in kNN
        S=s';
        [dim,N] = size(S);
        rrw = S(:);
        [kNNgraphlength, Graph] = kNNgraphmex(rrw, N, dim, kneighbors, 1);
        Graph = reshape(Graph, kneighbors+1, N );
end

if nargout>0
    N=size(s,1);
    G=sparse(N,N);
    for i = 1 : size(Graph,2)
    for j = 1 : size(Graph,1)   % k+1
         G(i,Graph(j,i))=1;
         G(Graph(j,i),i)=1;
    end
    end
    % G=0.5*(G+G');
    G=G-diag(diag(G));
end

if plotit
hold on
for i = 1 : size(Graph,2)
for j = 1 : size(Graph,1)
     line(s([i,Graph(j,i)],1),...
          s([i,Graph(j,i)],2),...
          s([i,Graph(j,i)],3),'Color','red');   
end
end
end

end
