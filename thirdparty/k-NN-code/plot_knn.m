function [Graph]=plot_knn(X, kneighbors)

% Plot k-NN graph (2D or 3D)
%
% Written by (c)Jose Costa (jcosta@umich.edu), version 1.0, 2004


plot_graph_option = 1;

[dim, N] = size(X);


rrw = X(:);
[kNNgraphlength, Graph] = kNNgraphmex(rrw, N, dim, kneighbors, 1);

Graph = reshape(Graph, kneighbors+1, N );

rr = X';


mark='o';
    lwidth=1;
    marksize=2.5;
    markfacecol='r';
    markedgecol='k';
    lcol='-r';
    
C={ mark, ...
    'markersize',marksize,...
    'markerfacecolor',markfacecol,...
    'markeredgecolor',markedgecol};
%C2={'-r','linewidth',lwidth};

if plot_graph_option==1 && dim==2
%    figure
    hold on
     for i = 1 : size(Graph,2)
         for j = 1 : size(Graph,1) 
             line(rr([i,Graph(j,i)],1),rr([i,Graph(j,i)],2))
         end
     end
    plot(rr(:,1),rr(:,2),'ro')
    %~ Uncomment following 3 lines if you want to see the point indices displayed
    %~ for i = 1 : N
    %~     text(rrw(i,1),rrw(i,2),int2str(i))
    %~ end
end
if plot_graph_option==1 && dim==3
    figure
    hold on
     for i = 1 : size(Graph,2)
         for j = 1 : size(Graph,1) 
             line(rr([i,Graph(j,i)],1),rr([i,Graph(j,i)],2), rr([i,Graph(j,i)],3))
         end
     end
    plot3(rr(:,1),rr(:,2),rr(:,3),C{:})
    % Uncomment following 3 lines if you want to see the point indices displayed
    %~ for i = 1 : N
    %~     text(rrw(i,1),rrw(i,2),rrw(i,3),int2str(i))
    %~ end
end

grid on, axis equal