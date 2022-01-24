function [Lineage, Ptime,Cell_dist] = Lineage_Ptime(W,No_cluster,cluster_label,root_cluster,root_cell,latent,reverse)
% Computing lineage and pseudotime.
%
% Input
%   -W:             Cel-to-cell similarity matrix
%   -No_cluster:    Number of clusters inferred by SoptSC
%   -cluster_label: Cluster labels from SoptSC
%   -root_cluster:  Set up root cluster; If root_cluster is set
%                   as 0, SoptSC automatically generates root cluster.
%   -root_cell:     Root cell specified by user; If root_cell is not
%                   provided (root_cell = 0), SoptSC will generate the 
%                   root cell in an unsupervised manner
%                   
%   -latent:        Latent space from SoptSC
%   -reverse:       Parameter to infer whether to reverse the lineage tree;
%                   Set reverse = 1 if you want to reverse the lineage
%                   tree. Otherwise set reverse = 0.
% Output
%   -Lineage:       Lineage identified by SoptSC.
%   -Ptime:         Pseudotime inferred by SoptSC.
%   -Cell_dist:     Cell-to-cell distance on graph.

%% Inferring psudotime
low_dis = squareform(pdist(latent));
No_cell = size(W,1);
W_graph1 = zeros(size(W));
W_graph1(W>0) = 1;

W_graph1(1:No_cell+1:end) = 0;

CC_Graph = graph(W_graph1);

% shortest path between cells
Short_pathd1 = distances(CC_Graph);

% initial root cell if root_cell == 0
if root_cell == 0
%     root_cell0 = 1;
    [aa1,aa2] = find(low_dis == max(low_dis(:)));
    root_cell0 = aa1(1);
else
    root_cell0 = root_cell;
end


%% lineage using binary graph
W_graph = W_graph1;

[nComponents,sizes,members] = networkComponents(W_graph);

AA = graph(W_graph);
Shortest_path = distances(AA);

Vertex = zeros(nComponents,2);
if nComponents >1
    for i = 1:nComponents
        memberi = members{i};
        ddi = low_dis(members{i},members{i});
        [ia,ib] = find(ddi==max(ddi(:)));
        Vertex(i,:) = [memberi(ia(1)) memberi(ib(1))];
    end
end


if nComponents >1
    root_zz = 1;
    for k1 = 1:nComponents
        if ismember(root_cell0,members{k1})
            root_zz = k1;
            break;
        end
    end
    root_zz1 = setdiff(1:nComponents,root_zz);
    
    % find closed point in root clique to root cell in Vertex(root_cell,:)
    ddroot = [Shortest_path(root_cell0,Vertex(root_zz,1)) Shortest_path(root_cell0,Vertex(root_zz,2))];
    [~,i0] = max(ddroot);
    start_cell = Vertex(root_zz,i0);
    uu1 = setdiff(Vertex(:),Vertex(root_zz,:));
    
    for ii = 1:nComponents-1
        if ii == 1
            vv = low_dis(start_cell,uu1);
            [ia,ib] = find(vv==min(vv(:)));
            W_graph(start_cell(ia),uu1(ib)) = 2;
            display([start_cell(ia),uu1(ib)]);
            
            for kki = 1:length(root_zz1)
                if ismember(uu1(ib),members{root_zz1(kki)})
                    start_cell = union(start_cell,Vertex(root_zz1(kki),:));
                    root_zz1 = setdiff(root_zz1,root_zz1(kki));
                    break;
                end
            end
            uu1 = setdiff(uu1,start_cell);
        else
            uu1 = setdiff(uu1,start_cell);
            vv = low_dis(start_cell,uu1);
            [ia,ib] = find(vv==min(vv(:)));
            W_graph(start_cell(ia),uu1(ib)) = 2;
            display([start_cell(ia),uu1(ib)]);
            
            for kki = 1:length(root_zz1)
                if ismember(uu1(ib),members{root_zz1(kki)})
                    start_cell = union(start_cell,Vertex(root_zz1(kki),:));
                    root_zz1 = setdiff(root_zz1,root_zz1(kki));
                    break;
                end
            end
        end
        display(root_zz1);
        display(start_cell);
        display(uu1);
    end
end


W_graph = (W_graph + W_graph')./2;



CC_Graph = graph(W_graph);
% figure;
% plot(CC_Graph,'Marker','s','MarkerSize',5,'Layout','force','NodeLabel',[],...
%     'NodeCData',true_labs,'NodeColor','flat','LineWidth',1);

% shortest path between cells
Short_pathd = distances(CC_Graph);

% Cluster to Cluster (C-to-C) graph (compute the distance)
CC_adjacent = zeros(No_cluster);
for i = 1:No_cluster
    for j = i+1:No_cluster
        sub_matrix = Short_pathd(cluster_label==i,cluster_label==j);
        [p,q] = size(sub_matrix);
        CC_adjacent(i,j) = sum(sub_matrix(:))./(p.*q);
    end
end

% construct CC graph
CC_Graph = graph(CC_adjacent,'upper');

if root_cluster ~= 0
    [Tree,pred] = minspantree(CC_Graph,'Root',root_cluster);
    % if root_cell = [], then infer root_cell for given root_cluster
    if root_cell == 0
        rootcc_idx = find(cluster_label==root_cluster);
        tau_score = zeros(1,length(rootcc_idx));
        Ave_ptime_cc = zeros(length(rootcc_idx),No_cluster);
        for jj = 1:length(rootcc_idx)
            [~,Ptimejj] = sort(Short_pathd1(rootcc_idx(jj),:));
            % Average ptime for each cluster
            for kk = 1:No_cluster
                Ave_ptime_cc(jj,kk) = mean(Ptimejj(cluster_label==kk));
            end
            Ave_ptime_cc(jj,:) = Ave_ptime_cc(jj,:)./max(Ave_ptime_cc(jj,:));
            tau_score(jj) = corr((pred'+1)./max(pred'+1),Ave_ptime_cc(jj,:)','type','Kendall');
        end
        [~,root_cell1] = max(tau_score);
        root_cell = rootcc_idx(root_cell1);
        display('Inferred root cell is:');
        display(root_cell);
%        plot(tau_score);
    end
    
else
    if root_cell ~= 0
        root_cluster = cluster_label(root_cell);
        [Tree,pred] = minspantree(CC_Graph,'Root',root_cluster);
    else
        % find the root cluster
        [a1,a2] = find(CC_adjacent == max(CC_adjacent(:)));
        display('Root cluster candidates:');
        display([a1 a2]);
        % root_cell = [], infer root_cell based on inferred root_cluster
        if reverse == 1
            [Tree,pred] = minspantree(CC_Graph,'Root',a1);
            
            if root_cell == 0
                rootcc_idx = find(cluster_label==a1);
                tau_score = zeros(1,length(rootcc_idx));
                Ave_ptime_cc = zeros(length(rootcc_idx),No_cluster);
                for jj = 1:length(rootcc_idx)
                    [~,Ptimejj] = sort(Short_pathd1(rootcc_idx(jj),:));
                    % Average ptime for each cluster
                    for kk = 1:No_cluster
                        Ave_ptime_cc(jj,kk) = mean(Ptimejj(cluster_label==kk));
                    end
                    Ave_ptime_cc(jj,:) = Ave_ptime_cc(jj,:)./max(Ave_ptime_cc(jj,:));
                    tau_score(jj) = corr((pred'+1)./max(pred'+1),Ave_ptime_cc(jj,:)','type','Kendall');
                end
                [~,root_cell1] = max(tau_score);
                root_cell = rootcc_idx(root_cell1);
                display('Inferred root cell is:');
                display(root_cell);
%                plot(tau_score);
            end
        else
            [Tree,pred] = minspantree(CC_Graph,'Root',a2);
            
            if root_cell == 0
                rootcc_idx = find(cluster_label==a2);
                tau_score = zeros(1,length(rootcc_idx));
                Ave_ptime_cc = zeros(length(rootcc_idx),No_cluster);
                for jj = 1:length(rootcc_idx)
                    [~,Ptimejj] = sort(Short_pathd1(rootcc_idx(jj),:));
                    % Average ptime for each cluster
                    for kk = 1:No_cluster
                        Ave_ptime_cc(jj,kk) = mean(Ptimejj(cluster_label==kk));
                    end
                    Ave_ptime_cc(jj,:) = Ave_ptime_cc(jj,:)./max(Ave_ptime_cc(jj,:));
                    tau_score(jj) = corr((pred'+1)./max(pred'+1),Ave_ptime_cc(jj,:)','type','Kendall');
                end
                [~,root_cell1] = max(tau_score);
                root_cell = rootcc_idx(root_cell1);
                display('Inferred root cell is:');
                display(root_cell);
%                plot(tau_score);
            end
            
            
        end
        
    end
end

% Compute CC MST
rootedTree = digraph(pred(pred~=0),find(pred~=0));
Lineage = pred;
figure;
plot(rootedTree);
set(gca,'xtick',[]);
set(gca,'ytick',[]);



%% plot ptime
Cell_dist = Short_pathd(root_cell,:);
[~,Ptime] = sort(Short_pathd1(root_cell,:));

% % plot pseudotime on the latent space
c = linspace(0,1,No_cell);
colormap parula;
cmap = colormap;
mymap = cmap(1:58,:);
colormap(mymap);
figure;
scatter(latent(Ptime,1),latent(Ptime,2),30,c,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);


cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
box on;

set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
