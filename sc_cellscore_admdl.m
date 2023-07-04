function [score]=sc_cellscore_admdl(X,genelist,tgsPos,tgsNeg,nbin,ctrl)
% Compute cell scores from a list of feature genes
%
% tgsPos - positive features (negative target marker genes)
% tgsNeg - negative features (negative target marker genes)
%
% see also: PKG.E_CELLSCORES, SC_CELLSCORE_UCELL 
%
% ref: AddModuleScore - https://github.com/satijalab/seurat/blob/master/R/utilities.R
% ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8271111/
% disp("Seurat provides a computationally efficient gene signature scoring function, 
% named AddModuleScore, originally proposed by Tirosh et al. [5]. ")

% The score is the average expression of a set of genes subtracted with 
% the average expression of a reference set of genes. The reference set is randomly sampled from the gene_pool for each binned expression value.
% This reproduces the approach in Seurat [Satija15] and has been implemented for Scanpy by Davide Cittaro.

if nargin<6, ctrl=5; end
if nargin<5, nbin=25; end
if nargin<4
    tgsNeg=[];
end
if nargin<3 || isempty(tgsPos)
    error('USAGE: >>[score]=sc_cellscore_admdl(X,genelist,tgsPos);')
    % tgsPos=["CD44","LY6C","KLRG1","CTLA","ICOS","LAG3"];
end

if ~any(matches(genelist, tgsPos,'IgnoreCase',true))
    score=NaN(size(X,2),1);
    warning('No feature genes found in GENELIST. NaN scores returned');
    return;
end

%genelist=upper(genelist);
%tgsPos=upper(tgsPos);
%tgsNeg=upper(tgsNeg);


X=sc_norm(X);
disp('Library-size normalization...done.')
X=log(X+1);
disp('Log(x+1) transformation...done.')



%idx=matches(genelist, tgsPos, 'IgnoreCase',true);
% if ~any(idx)
%     score=NaN(size(X,2),1);
%     warning('No feature genes found in GENELIST.');
%     return;
% end

[score]=i_calculate_score(X,genelist,tgsPos,1,nbin,ctrl);
if ~isempty(tgsNeg) && any(strlength(tgsNeg)>0)
    [s]=i_calculate_score(X,genelist,tgsNeg,-1,nbin,ctrl);
    score=score+s;    % add scores from negative markers (greater values for non-target cells)
end
end



function [score]=i_calculate_score(X,genelist,tgs,directtag,nbin,ctrl)
    if nargin<6, ctrl=5; end
    if nargin<5, nbin=25; end
    if nargin<4, directtag=1; end

    rng default
    % Initial stats
    cluster_lenght = size(X, 1);
    data_avg = mean(X, 2);
    [~, I] = sort(data_avg);

    % Sorting data
    data_avg = data_avg(I);
    gsorted = genelist(I);
    Xsorted = X(I, :);

    % Assigning bins by expression
    assigned_bin = diag(zeros(cluster_lenght));
    bin_size = cluster_lenght/nbin;
    for i = 1:nbin
        bin_match = data_avg <= data_avg(round(bin_size * i));
        pos_avail = (assigned_bin == 0);
        assigned_bin(pos_avail & bin_match) = i;
    end

    % Selecting bins of same expression
    idx=matches(gsorted, tgs,'IgnoreCase',true);
    selected_bins = unique(assigned_bin(idx));
    samebin_genes = gsorted(ismember(assigned_bin, selected_bins));
    ctrl_use = [];
    for i = 1:length(tgs)
        ctrl_use = [ctrl_use; ...
            randsample(samebin_genes, ctrl)];
    end
    ctrl_use = unique(ctrl_use);

    % Averaging expression
    ctrl_score = mean(Xsorted(matches(gsorted, ctrl_use,'IgnoreCase',true),:),1);
    features_score = mean(Xsorted(idx,:),1);

    % Scoring
    if directtag>0
        score = transpose(features_score - ctrl_score);
    else
        score = transpose(ctrl_score - features_score);
    end
end