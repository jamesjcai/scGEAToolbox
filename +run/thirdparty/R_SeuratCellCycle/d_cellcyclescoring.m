% Reading same input file
[X, genelist] = sc_readmtxfile('pbmc_small.mtx', 'genelist_pbmc_small.txt');

features = [
    "CD79B",
    "CD79A";
    "CD19";
    "CD180";
    "CD200";
    "CD3D";
    "CD2";
    "CD3E";
    "CD7";
    "CD8A";
    "CD14";
    "CD1C";
    "CD68";
    "CD9";
    "CD247"];

% Normalizing data by default in Seurat = log1p(libsize * 1e4)
X = log(((X./nansum(X)) * 1e4) + 1);

% Set seed
rng(1);

% Default arguments 
ctrl = 5;
nbin = 24;

% Initial stats
cluster_lenght = size(X, 1);
data_avg = mean(X, 2);
[~, sort_avg] = sort(data_avg);

% Sorting data
data_avg = data_avg(sort_avg);
genelist = genelist(sort_avg);
X = X(sort_avg, :);

% Assigning bins by expression
assigned_bin = diag(zeros(cluster_lenght));
bin_size = cluster_lenght/nbin;
for i = 1:nbin
    bin_match = data_avg <= data_avg(round(bin_size * i));
    pos_avail = (assigned_bin == 0);
    assigned_bin(pos_avail & bin_match) = i;
end

% Selecting bins of same expression
selected_bins = unique(assigned_bin(matches(genelist, features)));
samebin_genes = genelist(ismember(assigned_bin, selected_bins));
ctrl_use = [];
for i = 1:length(features)
    ctrl_use = [ctrl_use; randsample(samebin_genes, ctrl)];
end
ctrl_use = unique(ctrl_use);

% Averaging expression
ctrl_score = mean(X(matches(genelist, ctrl_use),:),1);
features_score = mean(X(matches(genelist, features),:),1)

% Scoring
score = features_score - ctrl_score;

csvwrite('out.csv', score)