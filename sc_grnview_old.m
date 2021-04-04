function sc_grnview_old(A,g,alpha)
if nargin<3, alpha=0.75; end
n=size(A,1);
if n>200, error('Smaller network is expected'); end
if nargin<2, g=string(1:n); end




B=e_transf(A,alpha);
G=digraph(B,g,'omitselfloops');
p = plot(G);
layout(p,'force');
G.Nodes.NodeColors = outdegree(G)-indegree(G);
p.NodeCData = G.Nodes.NodeColors;
G.Edges.LWidths = abs(7*G.Edges.Weight/max(G.Edges.Weight));
p.LineWidth = G.Edges.LWidths;
n=size(G.Edges,1);
cc=repmat([0 0.4470 0.7410],n,1);
cc(G.Edges.Weight<0,:)=repmat([0.8500, 0.3250, 0.0980],...
       sum(G.Edges.Weight<0),1);
p.EdgeColor=cc;
p.NodeFontSize=2*p.NodeFontSize;
%view(3)            
title('Single-cell Gene Regulatory Netowrk');
end


function A=e_transf(A,q)
% A - adjacency matrix
if nargin<2, q=0.95; end
dim=size(A);
if numel(dim)==2
    a=max(abs(A(:)));
    if a>0
        A=A./a;
        A=A.*(abs(A)>quantile(abs(A(:)),q));        
    end
elseif numel(dim)==3
    for k=1:dim(3)
        A(:,:,k)=e_transf(A(:,:,k),q);
    end
end
end