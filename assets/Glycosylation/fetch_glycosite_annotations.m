function fetch_glycosite_annotations()
%FETCH_GLYCOSITE_ANNOTATIONS  Rebuild glycosite_annotations.mat from UniProt.
%   Reproduces the data behind GLY.NODEANNOT from scratch: queries
%   UniProt's REST API for every gene in ASSETS/LIGAND_RECEPTOR's ligand and
%   receptor columns, parses each entry's CARBOHYD (glycosylation) and LIPID
%   (GPI-anchor) sequence features, and writes GLYCOSITE_ANNOTATIONS.MAT next
%   to this file.
%
%   Run this to REFRESH the asset against UniProt's current annotations (it
%   is updated continuously; the shipped asset is a dated snapshot - see its
%   source_note). Takes a few minutes: ~15 batched HTTP requests plus local
%   parsing of ~1400 entries. Requires internet access.
%
%   Not called by anything else in the toolbox; GLY.NODEANNOT only
%   reads the .mat file this produces.
%
% see also: GLY.NODEANNOT, GLY.LRWEIGHT

here = fileparts(mfilename('fullpath'));

% -- Gene list: every ligand/receptor symbol in the L-R database ----------
db = load(fullfile(here, '..', 'Ligand_Receptor', 'Ligand_Receptor.mat'), ...
    'ligand', 'receptor');
genes = unique(upper([string(db.ligand(:)); string(db.receptor(:))]));
fprintf('%d unique L-R gene symbols to query.\n', numel(genes));

% -- Query UniProt in batches of 90 genes (its OR-clause limit is 100) -----
opts = weboptions('Timeout', 60);
batchSize = 90;
nBatch = ceil(numel(genes)/batchSize);
raw = table();
for b = 1:nBatch
    lo = (b-1)*batchSize + 1;
    hi = min(b*batchSize, numel(genes));
    clause = strjoin("(gene:" + genes(lo:hi) + ")", " OR ");
    query = "(" + clause + ") AND organism_id:9606 AND reviewed:true";
    url = "https://rest.uniprot.org/uniprotkb/search";
    txt = webread(url, "query", query, ...
        "fields", "accession,gene_names,ft_carbohyd,ft_lipid", ...
        "format", "tsv", "size", 500, opts);
    tmpfile = [tempname, '.tsv'];
    fid = fopen(tmpfile, 'w');
    fwrite(fid, txt);
    fclose(fid);
    ropts = detectImportOptions(tmpfile, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    ropts.VariableNamingRule = 'preserve';
    part = readtable(tmpfile, ropts);
    delete(tmpfile);
    raw = [raw; part]; %#ok<AGROW>
    fprintf('  batch %d/%d: %d entries\n', b, nBatch, height(part));
end
writetable(raw, fullfile(here, 'uniprot_lr_raw.tsv'), 'FileType', 'text', 'Delimiter', '\t');

% -- Parse CARBOHYD / LIPID feature text into per-row site counts ---------
n = height(raw);
nNlink = zeros(n,1); nOGalNAc = zeros(n,1); nOFuc = zeros(n,1);
nOGlcNAc = zeros(n,1); nOOther = zeros(n,1); nGPI = zeros(n,1);
for i = 1:n
    t = raw.("Glycosylation")(i);
    if ~ismissing(t) && strlength(t) > 0
        notes = regexp(t, '/note="([^"]*)"', 'tokens');
        for k = 1:numel(notes)
            note = string(notes{k}{1});
            if contains(note, "N-linked", 'IgnoreCase', true)
                nNlink(i) = nNlink(i) + 1;
            elseif contains(note, "O-linked", 'IgnoreCase', true)
                if contains(note, "GalNAc", 'IgnoreCase', true)
                    nOGalNAc(i) = nOGalNAc(i) + 1;
                elseif contains(note, "Fuc", 'IgnoreCase', true)
                    nOFuc(i) = nOFuc(i) + 1;
                elseif contains(note, "GlcNAc", 'IgnoreCase', true)
                    nOGlcNAc(i) = nOGlcNAc(i) + 1;
                else
                    nOOther(i) = nOOther(i) + 1;
                end
            end
        end
    end
    tl = raw.("Lipidation")(i);
    if ~ismissing(tl) && strlength(tl) > 0
        lnotes = regexp(tl, '/note="([^"]*)"', 'tokens');
        for k = 1:numel(lnotes)
            if contains(lnotes{k}{1}, "GPI-anchor", 'IgnoreCase', true)
                nGPI(i) = 1;
            end
        end
    end
end

% -- Map each query gene to its UniProt row (first match wins) ------------
geneNameTokens = arrayfun(@(s) strsplit(s, ' '), raw.("Gene Names"), 'UniformOutput', false);
uniprotAcc = repmat("", numel(genes), 1);
found = false(numel(genes), 1);
n_nlinked = nan(numel(genes),1); n_o_galnac = nan(numel(genes),1);
n_o_fuc = nan(numel(genes),1); n_o_glcnac = nan(numel(genes),1);
n_o_other = nan(numel(genes),1); n_gpi = nan(numel(genes),1);
nMultiMatch = 0;
for gi = 1:numel(genes)
    hit = -1;
    for ri = 1:n
        if any(strcmpi(string(geneNameTokens{ri}), genes(gi)))
            if hit > 0, nMultiMatch = nMultiMatch + 1; else, hit = ri; end
        end
    end
    if hit > 0
        found(gi) = true;
        uniprotAcc(gi) = raw.Entry(hit);
        n_nlinked(gi) = nNlink(hit); n_o_galnac(gi) = nOGalNAc(hit);
        n_o_fuc(gi) = nOFuc(hit); n_o_glcnac(gi) = nOGlcNAc(hit);
        n_o_other(gi) = nOOther(hit); n_gpi(gi) = nGPI(hit);
    end
end
fprintf('matched %d of %d genes (%d had >1 UniProt hit, first kept)\n', ...
    sum(found), numel(genes), nMultiMatch);

T = table(genes, uniprotAcc, found, n_nlinked, n_o_galnac, n_o_fuc, n_o_glcnac, ...
    n_o_other, n_gpi, 'VariableNames', {'gene', 'uniprot_accession', 'found', ...
    'n_nlinked', 'n_o_galnac', 'n_o_fuc', 'n_o_glcnac', 'n_o_other', 'n_gpi'});
T.n_total = sum([T.n_nlinked, T.n_o_galnac, T.n_o_fuc, T.n_o_glcnac, T.n_o_other], ...
    2, 'omitnan');

source_note = "UniProt reviewed (Swiss-Prot) human proteome, CARBOHYD and LIPID " + ...
    "sequence features, fetched " + string(datetime('today', 'Format', 'yyyy-MM-dd')) + ...
    " via https://rest.uniprot.org/uniprotkb/search, restricted to gene symbols " + ...
    "appearing in assets/Ligand_Receptor/Ligand_Receptor.mat (n=" + numel(genes) + "). " + ...
    "Counts are documented site annotations, not predictions; a count of 0 means " + ...
    "no site is currently annotated in UniProt, not that the protein carries none. " + ...
    "Gene-symbol matching takes the first reviewed entry whose Gene Names field " + ...
    "contains the query symbol as any token; " + nMultiMatch + " of " + numel(genes) + ...
    " queries matched more than one entry and the first was kept.";

save(fullfile(here, 'glycosite_annotations.mat'), 'T', 'source_note');
fprintf('wrote %s\n', fullfile(here, 'glycosite_annotations.mat'));

end
