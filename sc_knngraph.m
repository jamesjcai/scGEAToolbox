function [Graph]=sc_knngraph(s,k, method)
if nargin<3, method=1; end
if nargin<2, k=4; end

switch method
    case 1
        [mIdx] = knnsearch(s,s,'K',k+1);
        Graph=mIdx';
    case 2
        pw1=fileparts(which(mfilename));
        pth=fullfile(pw1,'thirdparty/k-NN-code');
        addpath(pth);
        kneighbors=k; % number of neighbors in kNN
        S=s';
        [dim,N] = size(S);
        rrw = S(:);
        [kNNgraphlength, Graph] = kNNgraphmex(rrw, N, dim, kneighbors, 1);
        Graph = reshape(Graph, kneighbors+1, N );
end

hold on
for i = 1 : size(Graph,2)
 for j = 1 : size(Graph,1) 
     line(s([i,Graph(j,i)],1), s([i,Graph(j,i)],2), s([i,Graph(j,i)],3))
 end
end
end
