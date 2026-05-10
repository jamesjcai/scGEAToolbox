function net = sc_net_progeny(organism, options)
% SC_NET_PROGENY  Load the PROGENy 14-pathway signed gene-weight network.
%
% PROGENy (Pathway RespOnsive GENes) provides weighted gene sets for 14
% key signaling pathways derived from large-scale perturbation experiments.
% Unlike simple gene sets, weights encode the strength and direction of each
% gene's response to pathway activation.
%
% USAGE:
%   net = sc_net_progeny()                  % human, top 500 genes/pathway
%   net = sc_net_progeny("human")
%   net = sc_net_progeny("mouse")
%   net = sc_net_progeny("human", TopGenes=100)
%
% INPUT:
%   organism  - "human" (default) or "mouse"
%
% OUTPUT:
%   net - table with columns: source (pathway), target (gene), weight
%         Passes directly to sc_ulm().
%
% OPTIONS:
%   TopGenes (default 500) - top N genes per pathway ranked by |weight|
%
% PATHWAYS (14 total):
%   Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB,
%   p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT
%
% REFERENCE:
%   Schubert et al., Nature Communications 2018.
%   https://doi.org/10.1038/s41467-017-02391-6
%
% DATA SOURCE:
%   OmniPath web service (omnipathdb.org). Requires internet on first call;
%   result is cached locally in tempdir() for subsequent calls.
%   Alternatively, install decoupler-py and the Python client will be used.

arguments
    organism  (1,1) string = "human"
    options.TopGenes (1,1) double {mustBePositive} = 500
end

organism = lower(organism);
switch organism
    case {"human", "hs", "homo sapiens"}
        organism = "human";
        ncbi_id  = 9606;
    case {"mouse", "mm", "mus musculus"}
        organism = "mouse";
        ncbi_id  = 10090;
    otherwise
        error('sc_net_progeny:unknownOrganism', ...
            'organism must be "human" or "mouse". Got: %s', organism);
end

% --- Try cached file first ---
cacheFile = fullfile(tempdir(), sprintf('progeny_%s.mat', organism));
if isfile(cacheFile)
    s = load(cacheFile, 'net');
    net = s.net;
    net = i_filter_top(net, options.TopGenes);
    return
end

% --- Try Python decoupler if available ---
net = i_try_python_progeny(organism);
if ~isempty(net)
    save(cacheFile, 'net');
    net = i_filter_top(net, options.TopGenes);
    return
end

% --- Fetch from OmniPath REST API ---
fprintf('Downloading PROGENy network from OmniPath...\n');
url = sprintf([ ...
    'https://omnipathdb.org/queries/transcriptomics' ...
    '?resources=PROGENy&genesymbols=yes&organism=%d&format=tsv'], ncbi_id);

try
    opts = weboptions('Timeout', 30, 'ContentType', 'text');
    raw  = webread(url, opts);
    net  = i_parse_omnipath_tsv(raw, organism);
catch ME
    error('sc_net_progeny:downloadFailed', ...
        ['Could not download PROGENy data.\n' ...
         'Options:\n' ...
         '  1. Check internet connectivity.\n' ...
         '  2. Install decoupler-py: pip install decoupler\n' ...
         '  3. Provide a custom net table to sc_ulm().\n' ...
         'Error: %s'], ME.message);
end

save(cacheFile, 'net');
net = i_filter_top(net, options.TopGenes);
fprintf('PROGENy: %d pathways, %d gene-weight entries.\n', ...
    numel(unique(net.source)), height(net));

end


% -------------------------------------------------------------------------
function net = i_try_python_progeny(organism)
net = [];
try
    result = py.decoupler.op.progeny(organism=organism);
    source  = string(py.list(result{'source'}.tolist()));
    target  = string(py.list(result{'target'}.tolist()));
    weight  = double(py.array.array('d', result{'weight'}.tolist()));
    net = table(source(:), target(:), weight(:), ...
        'VariableNames', {'source', 'target', 'weight'});
catch
end
end


function net = i_parse_omnipath_tsv(raw, organism) %#ok<INUSD>
% Parse OmniPath TSV response into source/target/weight table.
lines  = strsplit(strtrim(raw), newline);
header = strsplit(lines{1}, '\t');
data   = cellfun(@(l) strsplit(l, '\t'), lines(2:end), 'UniformOutput', false);
data   = data(cellfun(@numel, data) == numel(header));

if isempty(data)
    error('sc_net_progeny:emptyResponse', ...
        'OmniPath returned no PROGENy entries. The API format may have changed.');
end

mat = vertcat(data{:});

% Locate relevant columns (OmniPath column names vary by endpoint version)
colSource  = i_find_col(header, {'source', 'pathway', 'set_name'});
colTarget  = i_find_col(header, {'target', 'genesymbol', 'gene_symbol'});
colWeight  = i_find_col(header, {'weight', 'mor', 'score', 'value'});

source = string(mat(:, colSource));
target = string(mat(:, colTarget));
weight = str2double(mat(:, colWeight));

valid = ~isnan(weight) & strlength(source) > 0 & strlength(target) > 0;
net = table(source(valid), target(valid), weight(valid), ...
    'VariableNames', {'source', 'target', 'weight'});
end


function col = i_find_col(header, candidates)
for c = candidates
    idx = find(strcmpi(header, c), 1);
    if ~isempty(idx)
        col = idx;
        return
    end
end
error('sc_net_progeny:columnNotFound', ...
    'Could not find a column matching: %s', strjoin(candidates, ', '));
end


function net = i_filter_top(net, topN)
if topN <= 0
    return
end
sources = unique(net.source);
keep = false(height(net), 1);
for k = 1:numel(sources)
    idx = net.source == sources(k);
    w   = abs(net.weight(idx));
    [~, order] = sort(w, 'descend');
    rows = find(idx);
    topRows = rows(order(1:min(topN, numel(order))));
    keep(topRows) = true;
end
net = net(keep, :);
end
