function proteins = sc_dock_gene2pdb(T_interactions, varargin)
% SC_DOCK_GENE2PDB  Map CellChat receptor gene names to PDB IDs.
%
% Bridges Module 2 (sc_dock_ccc) output to Module 3 (sc_dock_vina) input by
% converting receptor gene symbols (e.g. 'EGFR', 'VEGFR2') to PDB accession
% codes (e.g. '1IYT', '2GS6') suitable for structure download and docking.
%
% Lookup order:
%   1. Bundled receptor_reference.csv  (gene_name -> PDB_model column)
%   2. UniProt REST API fallback        (reviewed human entries, organism 9606)
%
% USAGE:
%   proteins = sc_dock_gene2pdb(ccc_result.T_interactions)
%   proteins = sc_dock_gene2pdb(ccc_result.T_interactions, 'top_n', 10)
%   proteins = sc_dock_gene2pdb(ccc_result.T_interactions, 'use_chain', true)
%
% INPUT:
%   T_interactions - table output from sc_dock_ccc (.T_interactions field),
%                   must contain a 'receptor' column of gene symbol strings.
%                   Rows should already be ranked by probability/delta_prob.
%
% OPTIONAL NAME-VALUE PAIRS:
%   'top_n'      - number of top-ranked interactions to consider (default 20)
%   'use_chain'  - if true, return full 'XXXX.C' chain identifiers instead of
%                  bare 4-character PDB IDs (default false)
%   'ref_file'   - path to an alternative gene->PDB reference CSV with columns
%                  'protein_name' and 'PDB_model' (default: bundled
%                  receptor_reference.csv co-located with this file)
%
% OUTPUT:
%   proteins - string array of unique PDB IDs (or PDB.chain strings when
%              use_chain=true), ready to pass to sc_dock_vina.
%              Genes with no PDB match are silently skipped with a warning.
%
% EXAMPLE:
%   ccc     = sc_dock_ccc(X_norm, g, c_cell_type);
%   proteins = sc_dock_gene2pdb(ccc.T_interactions, 'top_n', 15);
%   result  = sc_dock_vina(proteins, "fda");
%
% See also: SC_DOCK_CCC, SC_DOCK_VINA, SC_DOCK

p = inputParser;
addRequired(p,  'T_interactions', @istable);
addParameter(p, 'top_n',     20,    @(x) isscalar(x) && x >= 1);
addParameter(p, 'use_chain', false, @islogical);
addParameter(p, 'ref_file',  '',    @ischar);
parse(p, T_interactions, varargin{:});
opt = p.Results;

if ~ismember('receptor', T_interactions.Properties.VariableNames)
    error('sc_dock_gene2pdb:NoReceptorCol', ...
        'T_interactions must have a ''receptor'' column (output of sc_dock_ccc).');
end

top_n     = min(opt.top_n, height(T_interactions));
rec_genes = unique(T_interactions.receptor(1:top_n));

% Load reference table
ref_map = i_load_ref(opt.ref_file);

proteins  = strings(0);
no_match  = strings(0);

for i = 1:numel(rec_genes)
    gene = char(rec_genes(i));
    pdb_chain = '';

    if isKey(ref_map, gene)
        % CSV may list multiple models separated by "; " — take the first
        entries   = strsplit(ref_map(gene), '; ');
        pdb_chain = strtrim(entries{1});   % e.g. "1IYT.A"
    else
        % UniProt fallback — returns bare PDB ID, no chain
        pdb_chain = i_gene2pdb_uniprot(gene);
        if ~isempty(pdb_chain)
            fprintf('[sc_dock_gene2pdb] UniProt mapped %s -> %s\n', gene, pdb_chain);
        end
    end

    if isempty(pdb_chain)
        no_match(end+1) = string(gene); %#ok<AGROW>
        continue
    end

    if opt.use_chain
        % Keep full "XXXX.C" form
        proteins(end+1) = string(pdb_chain); %#ok<AGROW>
    else
        % Strip chain suffix: "1IYT.A" -> "1IYT"
        if contains(pdb_chain, '.')
            proteins(end+1) = string(extractBefore(pdb_chain, '.')); %#ok<AGROW>
        else
            proteins(end+1) = string(pdb_chain); %#ok<AGROW>
        end
    end
end

proteins = unique(proteins);

if ~isempty(no_match)
    warning('sc_dock_gene2pdb:NoMatch', ...
        'No PDB found for %d gene(s): %s', ...
        numel(no_match), strjoin(no_match, ', '));
end

fprintf('[sc_dock_gene2pdb] %d receptor gene(s) -> %d unique PDB ID(s).\n', ...
    numel(rec_genes), numel(proteins));
end

% =========================================================================
% Local helpers
% =========================================================================

function ref_map = i_load_ref(ref_file)
% Load gene->PDB mapping from receptor_reference.csv into a containers.Map.
ref_map = containers.Map('KeyType','char','ValueType','char');
if isempty(ref_file)
    candidates = {
        fullfile(fileparts(mfilename('fullpath')), 'receptor_reference.csv');
        fullfile(pwd, 'receptor_reference.csv');
    };
    for k = 1:numel(candidates)
        if isfile(candidates{k}), ref_file = candidates{k}; break; end
    end
end
if isempty(ref_file) || ~isfile(ref_file)
    warning('sc_dock_gene2pdb:NoRefFile', ...
        'receptor_reference.csv not found. All lookups will use UniProt API.');
    return
end
T = readtable(ref_file, 'TextType', 'string');
if ~all(ismember({'protein_name','PDB_model'}, T.Properties.VariableNames))
    error('sc_dock_gene2pdb:BadRefFile', ...
        'ref_file must have columns ''protein_name'' and ''PDB_model''.');
end
for k = 1:height(T)
    if ~ismissing(T.PDB_model(k)) && T.PDB_model(k) ~= ""
        ref_map(char(T.protein_name(k))) = char(T.PDB_model(k));
    end
end
fprintf('[sc_dock_gene2pdb] Loaded %d entries from %s\n', ref_map.Count, ref_file);
end

function pdb_id = i_gene2pdb_uniprot(gene)
% Query UniProt REST API for the best reviewed human PDB structure for a gene.
pdb_id = '';
try
    url = sprintf( ...
        'https://rest.uniprot.org/uniprotkb/search?query=gene:%s+organism_id:9606+reviewed:true&format=tsv&fields=accession,xref_pdb', ...
        urlencode(gene));
    txt  = webread(url);
    lines = strsplit(strtrim(txt), newline);
    if numel(lines) < 2, return; end
    cols = strsplit(lines{2}, char(9));
    if numel(cols) >= 2 && ~isempty(strtrim(cols{2}))
        pdbs   = strsplit(strtrim(cols{2}), ';');
        pdb_id = strtrim(pdbs{1});
    end
catch
end
end
