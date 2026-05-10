function net = sc_net_panglaodb(organism, options)
% SC_NET_PANGLAODB  Build a PanglaoDB cell-type marker network for sc_ulm().
%
% Loads canonical cell-type marker genes from the bundled PanglaoDB database
% and returns them as a network table (source = cell type, target = gene,
% weight = marker score). Pass the result directly to sc_ulm() to score each
% cell for enrichment of each cell-type's marker gene set.
%
% USAGE:
%   net = sc_net_panglaodb()            % human, default settings
%   net = sc_net_panglaodb("mouse")
%   net = sc_net_panglaodb("human", MinGeneScore=1.5)
%
% INPUT:
%   organism - "human" (default) or "mouse"
%
% OUTPUT:
%   net - table with columns: source (cell type), target (gene), weight
%         (log2-fold enrichment score from PanglaoDB). Passes to sc_ulm().
%
% OPTIONS:
%   MinGeneScore (default 0) - keep genes with weight >= this value (0 = all)
%
% EXAMPLE:
%   net = sc_net_panglaodb("human");
%   [scores, pvals, terms] = sc_ulm(sce, net, MinGenes=3);
%   T = sc_rankby_group(scores, terms, sce.c_cluster_id, TopN=3);
%
% REFERENCE:
%   Franzen et al., Database 2019.
%   https://doi.org/10.1093/database/baz046

arguments
    organism  (1,1) string = "human"
    options.MinGeneScore (1,1) double = 0
end

organism = lower(organism);
switch organism
    case {"human", "hs", "homo sapiens"}
        markerFile = "markerlist_hs_panglaodb.txt";
        weightFile = "markerweight_hs.txt";
    case {"mouse", "mm", "mus musculus"}
        markerFile = "markerlist_mm_panglaodb.txt";
        weightFile = "markerweight_mm.txt";
    otherwise
        error('sc_net_panglaodb:unknownOrganism', ...
            'organism must be "human" or "mouse". Got: %s', organism);
end

pw1 = fileparts(fileparts(fileparts(mfilename('fullpath'))));  % toolbox root (up from ml_decoupler/ → external/)
markerPath = fullfile(pw1, 'external', 'fun_alona_panglaodb', markerFile);
weightPath = fullfile(pw1, 'external', 'fun_alona_panglaodb', weightFile);

if ~isfile(markerPath)
    error('sc_net_panglaodb:fileNotFound', ...
        'PanglaoDB marker file not found:\n  %s', markerPath);
end

% --- Load gene weights (gene -> score) ---
geneWeights = containers.Map('KeyType', 'char', 'ValueType', 'double');
if isfile(weightPath)
    fid = fopen(weightPath, 'r');
    line = fgetl(fid);
    while ischar(line)
        parts = strsplit(strtrim(line), '\t');
        if numel(parts) >= 2
            gene = upper(strtrim(parts{1}));
            w    = str2double(parts{2});
            if ~isnan(w)
                geneWeights(gene) = w;
            end
        end
        line = fgetl(fid);
    end
    fclose(fid);
end

% --- Parse cell-type → gene list ---
sourceCells = {};
targetGenes = {};
weights     = [];

fid = fopen(markerPath, 'r');
line = fgetl(fid);
while ischar(line)
    line = strtrim(line);
    if isempty(line)
        line = fgetl(fid);
        continue
    end
    parts = strsplit(line, '\t');
    if numel(parts) < 2
        line = fgetl(fid);
        continue
    end
    cellType = strtrim(parts{1});
    geneStr  = strtrim(parts{2});
    genes    = upper(strtrim(strsplit(geneStr, ',')));
    genes    = genes(strlength(genes) > 0);

    for g = 1:numel(genes)
        gene = genes{g};
        if isKey(geneWeights, gene)
            w = geneWeights(gene);
        else
            w = 1;
        end
        if w >= options.MinGeneScore
            sourceCells{end+1, 1} = cellType; %#ok<AGROW>
            targetGenes{end+1, 1} = gene;      %#ok<AGROW>
            weights(end+1, 1)     = w;          %#ok<AGROW>
        end
    end
    line = fgetl(fid);
end
fclose(fid);

if isempty(sourceCells)
    error('sc_net_panglaodb:empty', ...
        'No entries loaded. Check MinGeneScore threshold (%.2f).', ...
        options.MinGeneScore);
end

net = table(string(sourceCells), string(targetGenes), weights, ...
    'VariableNames', {'source', 'target', 'weight'});

n_types = numel(unique(net.source));
fprintf('PanglaoDB (%s): %d cell types, %d marker entries.\n', ...
    organism, n_types, height(net));

end
