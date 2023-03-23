function [A]=sc_knngraph(s,k,plotit,method)
%Generate KNN group network from cell embeddings
%
% input: S - cell embedding coordinates
% output: A - sparse adjacency matrix 
%

if nargin<4, method=1; end
if nargin<3, plotit=false; end
if nargin<2 || isempty(k), k=4; end

switch method
    case 1
        [mIdx] = knnsearch(s,s,'K',k+1);
        Graph=mIdx';
    case 2
        pw1=fileparts(mfilename('fullpath'));
        pth=fullfile(pw1,'+run','thirdparty','k-NN-code');
        if ~(ismcc || isdeployed), addpath(pth); end
        kneighbors=k; % number of neighbors in kNN
        S=s';
        [dim,N] = size(S);
        rrw = S(:);
        [~, Graph] = kNNgraphmex(rrw, N, dim, kneighbors, 1);
        Graph = reshape(Graph, kneighbors+1, N );
end

if nargout>0 || plotit
    N=size(s,1);
    A=zeros(N,N);
    for i = 1 : size(Graph,2)
    for j = 1 : size(Graph,1)   % k+1
         A(i,Graph(j,i))=1;
         A(Graph(j,i),i)=1;
    end
    end
    % G=0.5*(G+G');
    A=A-diag(diag(A));
    A=sparse(A);
end

if plotit
    hold on
    for i = 1 : size(Graph,2)
        for j = 1 : size(Graph,1)
            % if i~=Graph(j,i)
            if A(i,Graph(j,i))>0
            if size(s,2)>=3
             line(s([i,Graph(j,i)],1),...
                  s([i,Graph(j,i)],2),...
                  s([i,Graph(j,i)],3),'Color','red');   
            else
             line(s([i,Graph(j,i)],1),...
                  s([i,Graph(j,i)],2),'Color','red');
            end
            end
        end
    end
    hold off
end

end