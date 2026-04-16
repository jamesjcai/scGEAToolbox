function result = sc_dock_vina(proteins, compounds, varargin)
% SC_DOCK_VINA  Module 3 of scDock: AutoDock Vina molecular docking.
%
% Runs AutoDock Vina to estimate binding affinities between a set of
% target proteins and small-molecule compounds. Returns ranked drug
% candidates by predicted binding affinity.
%
% USAGE:
%   result = sc_dock_vina(proteins, compounds)
%   result = sc_dock_vina(proteins, compounds, 'exhaustiveness', 8, ...
%                         'vina_exe', '/path/to/vina')
%
% INPUTS:
%   proteins  - struct array or string array of protein identifiers.
%               Each element specifies one receptor. Accepts either:
%               (a) string array of PDB IDs (fetched from RCSB if file absent)
%               (b) struct array with fields .id (string) and .pdbqt (file path)
%
%   compounds - struct array or string array of compound identifiers.
%               Accepts either:
%               (a) string array of PubChem CIDs or CAS numbers
%                   (SDF fetched from PubChem if file absent)
%               (b) struct array with fields .id (string) and .sdf (file path)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'outdir'          - output directory for docking files (default './dock_output')
%   'vina_exe'        - path to AutoDock Vina executable (default 'vina')
%   'obabel_exe'      - path to OpenBabel executable (default 'obabel')
%   'exhaustiveness'  - docking exhaustiveness 1-32 (default 8)
%   'n_modes'         - number of binding modes to generate (default 9)
%   'box_size'        - [x y z] docking box size in Angstroms
%                       (default [30 30 30]; global docking uses protein extent)
%   'global_docking'  - use whole-protein bounding box (default true)
%   'ph'              - pH for hydrogen assignment via obabel (default 7.4)
%
% OUTPUT:
%   result  - struct with fields:
%     .T_results  - table ranked by affinity (kcal/mol, lower = better):
%                   protein_id, compound_id, affinity_kcal_mol,
%                   pose_file, log_file
%     .outdir     - output directory path
%     .n_docked   - number of successfully completed docking runs
%
% REQUIREMENTS:
%   - AutoDock Vina (vina) must be on the system PATH or provided via 'vina_exe'
%   - OpenBabel (obabel) must be on the system PATH or provided via 'obabel_exe'
%
% See also: SC_DOCK_RNA, SC_DOCK_CCC, SC_DOCK

p = inputParser;
addRequired(p, 'proteins');
addRequired(p, 'compounds');
addParameter(p, 'outdir',          './dock_output', @ischar);
addParameter(p, 'vina_exe',        'vina',          @ischar);
addParameter(p, 'obabel_exe',      'obabel',        @ischar);
addParameter(p, 'exhaustiveness',  8,               @(x) isscalar(x) && x >= 1 && x <= 32);
addParameter(p, 'n_modes',         9,               @(x) isscalar(x) && x >= 1);
addParameter(p, 'box_size',        [30 30 30],      @(x) isvector(x) && numel(x) == 3);
addParameter(p, 'global_docking',  true,            @islogical);
addParameter(p, 'ph',              7.4,             @(x) isscalar(x) && x > 0);
parse(p, proteins, compounds, varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Expand 'fda' shorthand to CAS numbers from bundled fda.txt
% -------------------------------------------------------------------------
if isstring(compounds) && isscalar(compounds) && compounds == "fda"
    compounds = i_load_fda_compounds();
    fprintf('[sc_dock_vina] Loaded %d FDA compounds from fda.txt.\n', numel(compounds));
end

% -------------------------------------------------------------------------
% Setup output directory
% -------------------------------------------------------------------------
if ~exist(opt.outdir, 'dir')
    mkdir(opt.outdir);
end

% -------------------------------------------------------------------------
% Validate external tool availability
% -------------------------------------------------------------------------
i_check_exe(opt.vina_exe,   'AutoDock Vina');
i_check_exe(opt.obabel_exe, 'OpenBabel');

% -------------------------------------------------------------------------
% Resolve protein list to struct array with .id and .pdbqt fields
% -------------------------------------------------------------------------
proteins = i_resolve_proteins(proteins, opt.outdir, opt.obabel_exe, opt.ph);

% -------------------------------------------------------------------------
% Resolve compound list to struct array with .id and .sdf fields, then
% convert to PDBQT
% -------------------------------------------------------------------------
compounds = i_resolve_compounds(compounds, opt.outdir, opt.obabel_exe, opt.ph);

% -------------------------------------------------------------------------
% Run docking for each protein x compound pair
% -------------------------------------------------------------------------
n_p = numel(proteins);
n_c = numel(compounds);
n_total = n_p * n_c;
fprintf('[sc_dock_vina] Starting %d docking runs (%d proteins x %d compounds)...\n', ...
    n_total, n_p, n_c);

prot_ids  = strings(n_total, 1);
comp_ids  = strings(n_total, 1);
affinities = nan(n_total, 1);
pose_files = strings(n_total, 1);
log_files  = strings(n_total, 1);

run_idx = 0;
for pi = 1:n_p
    prot = proteins(pi);
    for ci = 1:n_c
        comp = compounds(ci);
        run_idx = run_idx + 1;
        fprintf('[sc_dock_vina]  (%d/%d) %s vs %s\n', ...
            run_idx, n_total, prot.id, comp.id);

        out_stem = fullfile(opt.outdir, ...
            sprintf('%s__%s', i_safe_name(prot.id), i_safe_name(comp.id)));
        pose_file = [out_stem '_out.pdbqt'];
        log_file  = [out_stem '_log.txt'];

        prot_ids(run_idx)  = prot.id;
        comp_ids(run_idx)  = comp.id;
        pose_files(run_idx) = string(pose_file);
        log_files(run_idx)  = string(log_file);

        % Compute docking box from protein structure
        if opt.global_docking
            [center, box_sz] = i_protein_bbox(prot.pdbqt);
        else
            center = [0 0 0];   % caller-specified center not yet exposed
            box_sz = opt.box_size;
        end

        % Build Vina command (--log removed in Vina v1.2+; capture stdout instead)
        cmd = sprintf(['"%s" --receptor "%s" --ligand "%s" ' ...
            '--center_x %.3f --center_y %.3f --center_z %.3f ' ...
            '--size_x %.1f --size_y %.1f --size_z %.1f ' ...
            '--exhaustiveness %d --num_modes %d ' ...
            '--out "%s"'], ...
            opt.vina_exe, prot.pdbqt, comp.pdbqt, ...
            center(1), center(2), center(3), ...
            box_sz(1), box_sz(2), box_sz(3), ...
            opt.exhaustiveness, opt.n_modes, ...
            pose_file);

        [status, cmdout] = system(cmd);
        % Write stdout to log file (replaces --log behaviour)
        fid = fopen(log_file, 'w');
        if fid ~= -1
            fprintf(fid, '%s', cmdout);
            fclose(fid);
        end
        if status ~= 0
            warning('sc_dock_vina:VinaFailed', ...
                'Vina failed for %s vs %s:\n%s', prot.id, comp.id, cmdout);
            continue
        end

        % Parse best affinity from log
        affinities(run_idx) = i_parse_vina_log(log_file);
    end
end

% -------------------------------------------------------------------------
% Build results table, rank by affinity (most negative = tightest binder)
% -------------------------------------------------------------------------
T_results = table(prot_ids, comp_ids, affinities, pose_files, log_files, ...
    'VariableNames', {'protein_id','compound_id','affinity_kcal_mol', ...
                      'pose_file','log_file'});
T_results = T_results(~isnan(T_results.affinity_kcal_mol), :);
T_results = sortrows(T_results, 'affinity_kcal_mol', 'ascend');

if height(T_results) > 0
    fprintf('[sc_dock_vina] Done. %d/%d runs succeeded. Best affinity: %.2f kcal/mol.\n', ...
        height(T_results), n_total, T_results.affinity_kcal_mol(1));
else
    warning('sc_dock_vina:NoDockingResults', ...
        'All %d docking runs failed. Check receptor PDBQT files and Vina installation.', n_total);
end

result.T_results  = T_results;
result.outdir     = opt.outdir;
result.n_docked   = height(T_results);
end

% =========================================================================
% Local helpers
% =========================================================================

function i_check_exe(exe, name)
[status, ~] = system(sprintf('"%s" --version', exe));
% Vina prints to stderr so also try --help
if status ~= 0
    [status2, ~] = system(sprintf('"%s" --help', exe));
    if status2 ~= 0
        error('sc_dock_vina:ExeNotFound', ...
            '%s not found at "%s". Add it to PATH or specify via parameter.', ...
            name, exe);
    end
end
end

function proteins = i_resolve_proteins(proteins, outdir, obabel_exe, ph)
if isstring(proteins) || iscellstr(proteins)
    proteins = string(proteins);
    tmp = struct('id', cell(numel(proteins),1), 'pdbqt', cell(numel(proteins),1));
    for i = 1:numel(proteins)
        id_str = char(proteins(i));
        tmp(i).id = id_str;
        % Parse optional chain suffix: "1A22.B" -> pdb_base="1A22", chain="B"
        if contains(id_str, '.')
            parts    = strsplit(id_str, '.');
            pdb_base = upper(parts{1});
            chain    = upper(parts{2});
            safe_id  = [pdb_base '_' chain];
        else
            pdb_base = upper(id_str);
            chain    = '';
            safe_id  = pdb_base;
        end
        pdb_file   = fullfile(outdir, [safe_id '.pdb']);
        pdbqt_file = fullfile(outdir, [safe_id '.pdbqt']);
        if ~isfile(pdb_file)
            i_fetch_pdb(pdb_base, pdb_file, chain);
        end
        if ~isfile(pdbqt_file)
            i_pdb2pdbqt(pdb_file, pdbqt_file, obabel_exe, ph);
        end
        tmp(i).pdbqt = pdbqt_file;
    end
    proteins = tmp;
else
    % Already a struct array; ensure pdbqt files exist
    for i = 1:numel(proteins)
        if ~isfield(proteins(i), 'pdbqt') || ~isfile(proteins(i).pdbqt)
            pdbqt_file = fullfile(outdir, [i_safe_name(proteins(i).id) '.pdbqt']);
            if isfield(proteins(i), 'pdb') && isfile(proteins(i).pdb)
                i_pdb2pdbqt(proteins(i).pdb, pdbqt_file, obabel_exe, ph);
            else
                error('sc_dock_vina:MissingProtein', ...
                    'Cannot find/prepare PDBQT for protein %s.', proteins(i).id);
            end
            proteins(i).pdbqt = pdbqt_file;
        end
    end
end
end

function compounds = i_resolve_compounds(compounds, outdir, obabel_exe, ph)
if isstring(compounds) || iscellstr(compounds)
    compounds = string(compounds);
    tmp = struct('id', cell(numel(compounds),1), 'pdbqt', cell(numel(compounds),1));
    for i = 1:numel(compounds)
        tmp(i).id   = char(compounds(i));
        sdf_file    = fullfile(outdir, [i_safe_name(char(compounds(i))) '.sdf']);
        pdbqt_file  = fullfile(outdir, [i_safe_name(char(compounds(i))) '.pdbqt']);
        if ~isfile(sdf_file)
            i_fetch_pubchem(compounds(i), sdf_file);
        end
        if ~isfile(pdbqt_file)
            i_sdf2pdbqt(sdf_file, pdbqt_file, obabel_exe, ph);
        end
        tmp(i).pdbqt = pdbqt_file;
    end
    compounds = tmp;
else
    for i = 1:numel(compounds)
        if ~isfield(compounds(i), 'pdbqt') || ~isfile(compounds(i).pdbqt)
            pdbqt_file = fullfile(outdir, [i_safe_name(compounds(i).id) '.pdbqt']);
            if isfield(compounds(i), 'sdf') && isfile(compounds(i).sdf)
                i_sdf2pdbqt(compounds(i).sdf, pdbqt_file, obabel_exe, ph);
            else
                error('sc_dock_vina:MissingCompound', ...
                    'Cannot find/prepare PDBQT for compound %s.', compounds(i).id);
            end
            compounds(i).pdbqt = pdbqt_file;
        end
    end
end
end

function i_fetch_pdb(pdb_id, out_file, chain_id)
% Download PDB from RCSB. If chain_id is non-empty, keep only that chain's
% ATOM/TER records (replicates download_pdb_chains_from_csv.py chain extraction).
if nargin < 3, chain_id = ''; end
base_id = upper(char(pdb_id));
url = sprintf('https://files.rcsb.org/download/%s.pdb', base_id);
fprintf('    Downloading PDB %s...\n', base_id);
tmp_file = [out_file '.download.tmp'];
try
    websave(tmp_file, url);
catch ME
    error('sc_dock_vina:PDBFetch', ...
        'Failed to download PDB %s: %s', base_id, ME.message);
end
if isempty(chain_id)
    movefile(tmp_file, out_file);
    return
end
% Filter to requested chain (PDB format: chain ID is column 22)
fid_in  = fopen(tmp_file,  'r');
fid_out = fopen(out_file, 'w');
n_kept = 0;
while ~feof(fid_in)
    line = fgetl(fid_in);
    if ~ischar(line), break; end
    rec = strtrim(line(1:min(6, end)));
    if (strcmp(rec, 'ATOM') || strcmp(rec, 'HETATM')) && ...
            numel(line) >= 22 && upper(line(22)) == upper(chain_id(1))
        fprintf(fid_out, '%s\n', line);
        n_kept = n_kept + 1;
    elseif strcmp(rec, 'TER') || strcmp(rec, 'END')
        fprintf(fid_out, '%s\n', line);
    end
end
fclose(fid_in);
fclose(fid_out);
delete(tmp_file);
fprintf('    Extracted chain %s: %d ATOM records -> %s\n', chain_id, n_kept, out_file);
end

function i_fetch_pubchem(compound_id, out_file)
% Try PubChem CID first, then name search
cid = str2double(char(compound_id));
if ~isnan(cid)
    url = sprintf( ...
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%d/SDF', cid);
else
    url = sprintf( ...
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/SDF', ...
        urlencode(char(compound_id)));
end
fprintf('    Downloading SDF for %s from PubChem...\n', compound_id);
try
    websave(out_file, url);
catch ME
    error('sc_dock_vina:PubChemFetch', ...
        'Failed to download compound %s: %s', compound_id, ME.message);
end
end

function i_pdb2pdbqt(pdb_file, pdbqt_file, obabel_exe, ph)
% Remove water/heteroatoms, add hydrogens at given pH, output rigid receptor PDBQT.
% -xr suppresses ROOT/ENDROOT/BRANCH/TORSDOF tags that Vina rejects in receptors.
cmd = sprintf('"%s" "%s" -O "%s" --addpolarh -p %.1f --partialcharge gasteiger -xr 2>&1', ...
    obabel_exe, pdb_file, pdbqt_file, ph);
[status, out] = system(cmd);
if status ~= 0
    warning('sc_dock_vina:ObabelFailed', ...
        'obabel protein conversion failed:\n%s', out);
end
end

function i_sdf2pdbqt(sdf_file, pdbqt_file, obabel_exe, ph)
cmd = sprintf('"%s" "%s" -O "%s" --addpolarh -p %.1f --partialcharge gasteiger --gen3d 2>&1', ...
    obabel_exe, sdf_file, pdbqt_file, ph);
[status, out] = system(cmd);
if status ~= 0
    warning('sc_dock_vina:ObabelFailed', ...
        'obabel compound conversion failed:\n%s', out);
end
end

function [center, box_sz] = i_protein_bbox(pdbqt_file)
% Read ATOM/HETATM records and compute bounding box + 10 Å padding
fid = fopen(pdbqt_file, 'r');
xyz = [];
while ~feof(fid)
    line = fgetl(fid);
    if length(line) >= 54 && (startsWith(line,'ATOM') || startsWith(line,'HETATM'))
        x = str2double(line(31:38));
        y = str2double(line(39:46));
        z = str2double(line(47:54));
        xyz(end+1, :) = [x y z]; %#ok<AGROW>
    end
end
fclose(fid);
if isempty(xyz)
    center = [0 0 0];
    box_sz = [30 30 30];
    return
end
mn = min(xyz, [], 1);
mx = max(xyz, [], 1);
center = (mn + mx) / 2;
box_sz = (mx - mn) + 20;  % 10 Å padding on each side
box_sz = max(box_sz, [20 20 20]);
end

function aff = i_parse_vina_log(log_file)
% Parse best (first) mode affinity from Vina log
aff = NaN;
if ~isfile(log_file), return; end
fid = fopen(log_file, 'r');
in_table = false;
while ~feof(fid)
    line = fgetl(fid);
    if contains(line, '-----+------------+----------+----------')
        in_table = true;
        continue
    end
    if in_table && ~isempty(strtrim(line))
        tokens = strsplit(strtrim(line));
        if numel(tokens) >= 2
            val = str2double(tokens{2});
            if ~isnan(val)
                aff = val;
                break
            end
        end
    end
end
fclose(fid);
end

function s = i_safe_name(s)
% Replace characters unsafe for filenames
s = regexprep(char(s), '[^a-zA-Z0-9_\-]', '_');
end

function cas_list = i_load_fda_compounds()
% Load CAS numbers from bundled fda.txt (one CAS per line).
candidates = {
    fullfile(fileparts(mfilename('fullpath')), 'fda.txt');
    fullfile(pwd, 'fda.txt');
};
f = '';
for k = 1:numel(candidates)
    if isfile(candidates{k}), f = candidates{k}; break; end
end
if isempty(f)
    error('sc_dock_vina:NoFDA', ...
        ['fda.txt not found. Download it from ' ...
         'https://github.com/Andrewneteye4343/scDock/tree/main/FDA_parent_compounds_202509 ' ...
         'and place it in the same folder as sc_dock_vina.m.']);
end
lines = strtrim(splitlines(fileread(f)));
cas_list = string(lines(~cellfun(@isempty, lines)));
end

function i_annotate_drug_info(cas_csv, drug_tsv)
% Annotate drug info table by joining CAS/CID CSV with a drug database TSV
% on the shared 'PubChem CID' column. Overwrites cas_csv in place.
%
% Replicates annotate_drug_info.py using native MATLAB table operations.
%
% INPUTS:
%   cas_csv  - path to CSV with at least columns 'CAS number' and 'PubChem CID'
%              (produced by i_fetch_pubchem / download_cas_pubchem.py)
%   drug_tsv - path to tab-delimited drug database with a 'PubChem CID' column
%              (e.g. DrugBank, ChEMBL export)
if ~isfile(cas_csv)
    error('sc_dock_vina:AnnotateDrug', 'cas_csv not found: %s', cas_csv);
end
if ~isfile(drug_tsv)
    error('sc_dock_vina:AnnotateDrug', 'drug_tsv not found: %s', drug_tsv);
end
fprintf('[sc_dock_vina] Loading %s\n', cas_csv);
T_cas  = readtable(cas_csv,  'TextType', 'string');
fprintf('[sc_dock_vina] Loading %s\n', drug_tsv);
T_drug = readtable(drug_tsv, 'FileType', 'text', 'Delimiter', '\t', ...
    'TextType', 'string');
% Normalise PubChem CID: strip trailing ".0", take first token before ";"
T_cas.("PubChem_CID")  = i_clean_cid(T_cas.("PubChem_CID"));
T_drug.("PubChem_CID") = i_clean_cid(T_drug.("PubChem_CID"));
% Left join on PubChem CID
T_out = outerjoin(T_cas, T_drug, 'Keys', 'PubChem_CID', ...
    'MergeKeys', true, 'Type', 'left');
% Put CAS number and PubChem CID first
vars     = T_out.Properties.VariableNames;
priority = {'CAS_number', 'PubChem_CID'};
other    = vars(~ismember(vars, priority));
T_out    = T_out(:, [priority(ismember(priority, vars)), other]);
n_matched = sum(ismember(T_out.("PubChem_CID"), T_drug.("PubChem_CID")));
fprintf('[sc_dock_vina] Matched %d / %d rows. Writing %s\n', ...
    n_matched, height(T_out), cas_csv);
writetable(T_out, cas_csv);
end

function cid = i_clean_cid(cid)
% Normalise PubChem CID strings: take first token before ";", strip ".0" suffix.
for k = 1:numel(cid)
    s = strtrim(char(cid(k)));
    if contains(s, ';'), s = strtrim(extractBefore(s, ';')); end
    if length(s) > 2 && strcmp(s(end-1:end), '.0'), s = s(1:end-2); end
    cid(k) = string(s);
end
end
