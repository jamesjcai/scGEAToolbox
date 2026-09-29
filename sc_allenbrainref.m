function ref = sc_allenbrainref(options)
%SC_ALLENBRAINREF Mouse whole-brain cell-type reference from the Allen Brain Cell Atlas.
%
%   ref = SC_ALLENBRAINREF() returns the Allen Brain Cell Atlas whole mouse
%   brain taxonomy (Yao et al., Nature 2023; about 4 million cells, 5,322
%   clusters) as a reference for SC_CELLTYPEANNOREF and
%   SC_ANNOTATECELLS(Method="reference"), at the level of its 338
%   subclasses - the level SAHA preloads.
%
%   ref = SC_ALLENBRAINREF(Level="class") gives the 34 classes instead.
%
%   DOWNLOAD ON FIRST USE. The atlas is not shipped with the toolbox: it is
%   CC BY-NC 4.0 (non-commercial use, with attribution to the Allen
%   Institute), and the toolbox is MIT. The first call downloads about
%   1.38 GB from the Allen Institute's public bucket, builds both levels,
%   caches them (about 30 MB) in CacheDir, and deletes the download. Later
%   calls load the cache. By using the reference you accept the Allen
%   Institute's terms; ref.License and ref.Citation say what they are.
%
%   Files downloaded, from s3://allen-brain-cell-atlas (no login):
%     mapmycells/WMB-taxonomy/20240831/precomputed_stats_ABC_revision_230821.h5
%         per-cluster mean expression, log2(CPM+1), 32,285 genes
%     metadata/WMB-taxonomy/20231215/views/cluster_annotation_term_with_counts.csv
%         the taxonomy and each cluster's cell count
%     metadata/WMB-10X/20231215/gene.csv
%         Ensembl ID to gene symbol
%
%   HOW IT IS BUILT. A subclass (or class) profile is the cell-weighted mean
%   of its clusters' mean profiles, so it equals the mean over its cells.
%   Genes are kept if some type expresses them at log2(CPM+1) >= 1.
%   Markers come from the means alone, because the published statistics
%   carry no detection rates: a gene is a marker of a type if the type's
%   mean is >= 1 and exceeds the cell-weighted mean of all other types by
%   >= 1 (two-fold), strongest 100 first. That differs from the Wilcoxon and
%   detection-rate filter SC_BUILDCELLTYPEREF uses on cell-level data, so
%   the marker test is somewhat less comparable here than the correlation.
%
%   Name-value options:
%     Level          "subclass" (default) or "class"
%     CacheDir       where the built reference is kept
%                    (prefdir/scgeatoolbox_addons/AllenBrainCellAtlas)
%     SourceDir      a folder already holding the three files; skips the
%                    download and leaves the files in place
%     AllowDownload  false errors instead of downloading when nothing is
%                    cached (true)
%     Rebuild        rebuild even when a cache exists (false)
%
%   Example:
%     ref = sc_allenbrainref();
%     sce = sc_annotatecells(sce, Method="reference", Reference=ref);
%
%   See also SC_CELLTYPEANNOREF, SC_BUILDCELLTYPEREF, SC_ANNOTATECELLS.

arguments
    options.Level (1, 1) string {mustBeMember(options.Level, ["subclass", "class"])} = "subclass"
    options.CacheDir (1, 1) string = fullfile(prefdir, 'scgeatoolbox_addons', 'AllenBrainCellAtlas')
    options.SourceDir (1, 1) string = ""
    options.AllowDownload (1, 1) logical = true
    options.Rebuild (1, 1) logical = false
end

cacheFile = in_cachefile(options.CacheDir, options.Level);
if isfile(cacheFile) && ~options.Rebuild
    S = load(cacheFile, 'ref');
    ref = S.ref;
    return
end

files = in_files();
if strlength(options.SourceDir) > 0
    srcDir = options.SourceDir;
    isTemporary = false;
else
    if ~options.AllowDownload
        error('sc_allenbrainref:NotCached', ...
            ['The Allen Brain Cell Atlas reference is not cached in %s. ', ...
             'Call sc_allenbrainref with AllowDownload=true to download it ', ...
             '(about 1.38 GB, once).'], options.CacheDir);
    end
    srcDir = fullfile(options.CacheDir, 'download');
    in_download(files, srcDir);
    isTemporary = true;
end

refs = in_build(srcDir, files);
if ~isfolder(options.CacheDir)
    mkdir(options.CacheDir);
end
for level = ["subclass", "class"]
    ref = refs.(level);
    save(in_cachefile(options.CacheDir, level), 'ref', '-v7');
end
if isTemporary
    % The built caches are what is kept; the 1.38 GB download is not.
    rmdir(srcDir, 's');
end
ref = refs.(options.Level);
end


function f = in_cachefile(cacheDir, level)
f = fullfile(cacheDir, "abc_wmb_ccn20230722_" + level + ".mat");
end


function files = in_files()
base = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/";
files = struct( ...
    "Stats", struct("Name", "precomputed_stats_ABC_revision_230821.h5", ...
        "Url", base + "mapmycells/WMB-taxonomy/20240831/precomputed_stats_ABC_revision_230821.h5"), ...
    "Terms", struct("Name", "cluster_annotation_term_with_counts.csv", ...
        "Url", base + "metadata/WMB-taxonomy/20231215/views/cluster_annotation_term_with_counts.csv"), ...
    "Genes", struct("Name", "gene.csv", ...
        "Url", base + "metadata/WMB-10X/20231215/gene.csv"));
end


function in_download(files, srcDir)
if ~isfolder(srcDir)
    mkdir(srcDir);
end
opts = weboptions('Timeout', 600);
for key = ["Terms", "Genes", "Stats"]
    target = fullfile(srcDir, files.(key).Name);
    if isfile(target)
        continue
    end
    fprintf('[allen] downloading %s...\n', files.(key).Name);
    partial = target + ".part";
    websave(partial, files.(key).Url, opts);
    movefile(partial, target);
end
end


function refs = in_build(srcDir, files)
statsFile = fullfile(srcDir, files.Stats.Name);
terms = readtable(fullfile(srcDir, files.Terms.Name), 'TextType', 'string', ...
    'Delimiter', ',');
genes = readtable(fullfile(srcDir, files.Genes.Name), 'TextType', 'string', ...
    'Delimiter', ',');

% The taxonomy: each cluster's supertype, subclass and class, by walking
% parent_term_label upward.
isCluster = terms.cluster_annotation_term_set_name == "cluster";
clusters = terms(isCluster, :);
[~, loc] = ismember(clusters.parent_term_label, terms.label);          % supertype
[~, loc] = ismember(terms.parent_term_label(loc), terms.label);        % subclass
subclassOf = terms.name(loc);
[~, loc2] = ismember(terms.parent_term_label(loc), terms.label);       % class
classOf = terms.name(loc2);
cellsPerCluster = double(clusters.number_of_cells);

% Rows of the statistics are clusters in the order cluster_to_row gives.
rowOf = jsondecode(char(h5read(statsFile, '/cluster_to_row')));
rows = zeros(height(clusters), 1);
for i = 1:height(clusters)
    rows(i) = rowOf.(clusters.label(i)) + 1;
end
ensembl = string(jsondecode(char(h5read(statsFile, '/col_names'))));
numGenes = numel(ensembl);

levels = struct("subclass", subclassOf, "class", classOf);
sums = struct();
for level = ["subclass", "class"]
    [~, ~, idx] = unique(levels.(level));
    sums.(level) = zeros(numGenes, max(idx));
end

% Accumulate cell-weighted sums of cluster means, a block of clusters at a
% time: the whole matrix is 1.4 GB in double.
blockSize = 500;
numRows = max(rows);
for first = 1:blockSize:numRows
    count = min(blockSize, numRows - first + 1);
    block = h5read(statsFile, '/sum', [1, first], [numGenes, count]);
    inBlock = find(rows >= first & rows < first + count);
    for level = ["subclass", "class"]
        [~, ~, idx] = unique(levels.(level));
        W = sparse(rows(inBlock) - first + 1, idx(inBlock), ...
            cellsPerCluster(inBlock), count, size(sums.(level), 2));
        sums.(level) = sums.(level) + block*W;
    end
end

[symbols, keepGene] = in_symbols(ensembl, genes);
refs = struct();
for level = ["subclass", "class"]
    [types, ~, idx] = unique(levels.(level));
    numCells = accumarray(idx, cellsPerCluster)';
    profile = sums.(level)(keepGene, :)./numCells;
    isExpressed = max(profile, [], 2) >= 1;
    profile = profile(isExpressed, :);
    geneNames = symbols(isExpressed);
    refs.(level) = struct( ...
        "Genes", geneNames, ...
        "Types", types(:)', ...
        "Profile", single(profile), ...
        "Markers", {in_markers(profile, numCells, geneNames)}, ...
        "NumCells", numCells, ...
        "Species", "mouse", ...
        "Source", "Allen Brain Cell Atlas, whole mouse brain CCN20230722, " + level, ...
        "License", "CC BY-NC 4.0 (non-commercial use, attribution to the Allen Institute)", ...
        "Citation", "Yao Z et al. A high-resolution transcriptomic and spatial atlas " + ...
            "of cell types in the whole mouse brain. Nature 624, 317-332 (2023). " + ...
            "doi:10.1038/s41586-023-06812-z");
end
end


function [symbols, keep] = in_symbols(ensembl, genes)
% Gene symbols for the Ensembl IDs; an ID without a symbol, and every
% repeat of a symbol after the first, is dropped.
[found, loc] = ismember(ensembl, genes.gene_identifier);
symbols = strings(numel(ensembl), 1);
symbols(found) = genes.gene_symbol(loc(found));
keep = found & strlength(symbols) > 0;
[~, firstOf] = unique(symbols, 'stable');
isFirst = false(numel(symbols), 1);
isFirst(firstOf) = true;
keep = keep & isFirst;
symbols = symbols(keep);
end


function markers = in_markers(profile, numCells, geneNames)
minMean = 1;
minDifference = 1;
numMarkers = 100;
numTypes = size(profile, 2);
total = profile*numCells(:);
markers = cell(1, numTypes);
for t = 1:numTypes
    restMean = (total - profile(:, t)*numCells(t))/(sum(numCells) - numCells(t));
    difference = profile(:, t) - restMean;
    candidate = find(profile(:, t) >= minMean & difference >= minDifference);
    [~, order] = sort(difference(candidate), 'descend');
    markers{t} = geneNames(candidate(order(1:min(numMarkers, numel(order)))));
end
end
