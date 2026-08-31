function [block, celltypes] = i_cellcontext(tissue, options)
%I_CELLCONTEXT Reference marker dictionary for one tissue, as prompt text.
%
%   [block, celltypes] = LLM.I_CELLCONTEXT(tissue)
%   [...] = LLM.I_CELLCONTEXT(tissue, Refresh=true, Verbose=true)
%
%   BLOCK is every cell type catalogued for TISSUE, one per line, formatted
%   as the CellTypeAI project formats it:
%
%       - Naive B cells: [CD19, MS4A1, CD79A, ...]
%
%   It is meant to be pasted into an annotation prompt under the heading
%   "CELL MARKERS FOR CONTEXT (Reference Dictionary)", giving the model a
%   fixed vocabulary to match a cluster's markers against instead of naming
%   a cell type from memory.
%
%   LLM.I_CELLTISSUES lists the tissues the dictionary covers, without
%   needing a network connection, so a chooser can be shown before the
%   markers are fetched.
%
%   BLOCK is "" when TISSUE is unset, unknown, or the dictionary cannot be
%   reached. That is not an error - the caller is expected to fall back to
%   a prompt without the reference section.
%
%   The dictionary is downloaded from the CellTypeAI project and cached
%   under prefdir. It is NOT shipped with this toolbox: CellTypeAI is
%   GPL-2.0 and scGEAToolbox is MIT, so redistributing it inside the
%   toolbox would put the two licences in conflict. Fetching it onto the
%   user's own machine at run time keeps them separate.
%
%   Source: https://github.com/rhdaw/CellTypeAI (GPL-2.0)
%
%   See also LLM.I_CELLTISSUES, LLM.E_CELLTYPEANNO, SC_ANNOTATECELLS.

arguments
    tissue (1,1) string = ""
    options.Refresh (1,1) logical = false
    options.Verbose (1,1) logical = false
end

block = "";
celltypes = strings(0, 1);

if strlength(strtrim(tissue)) == 0, return; end

recs = i_dictionary(options.Refresh, options.Verbose);
if isempty(recs), return; end

keep = strcmpi(string({recs.Tissue}), strtrim(tissue));
if ~any(keep), return; end

recs = recs(keep);
celltypes = strings(numel(recs), 1);
lines = strings(numel(recs), 1);
for k = 1:numel(recs)
    celltypes(k) = i_clean(recs(k).CellType);
    markers = cellfun(@i_clean, recs(k).Markers(:)');
    markers = markers(strlength(markers) > 0);
    lines(k) = "- " + celltypes(k) + ": [" + strjoin(markers, ", ") + "]";
end

drop = strlength(celltypes) == 0;
celltypes(drop) = [];
lines(drop) = [];
block = strjoin(lines, newline);
end



function recs = i_dictionary(refresh, verbose)
% The decoded dictionary, cached in memory and on disk. Returns [] rather
% than erroring when it cannot be obtained, so annotation degrades to a
% prompt without the reference section instead of failing outright.
persistent cached
if ~isempty(cached) && ~refresh
    recs = cached;
    return
end

url = "https://raw.githubusercontent.com/rhdaw/CellTypeAI/main/src/celltypeai/cell_context.json";
cachedir = fullfile(prefdir, 'scgeatoolbox_addons', 'celltypeai');
cachefile = fullfile(cachedir, 'cell_context.json');

if refresh || ~isfile(cachefile)
    if ~isfolder(cachedir), mkdir(cachedir); end
    try
        if verbose
            fprintf('[celltypeanno] downloading marker dictionary...\n');
        end
        websave(cachefile, url, weboptions('Timeout', 30));
    catch ME
        warning('llm:i_cellcontext:DownloadFailed', ...
            ['Could not download the cell-type marker dictionary (%s). ', ...
            'Annotation will run without the reference list.'], ME.message);
        recs = [];
        return
    end
end

try
    recs = i_validate(jsondecode(fileread(cachefile)));
catch ME
    % A truncated or malformed cache is worth deleting: the next call then
    % re-downloads instead of failing the same way forever.
    if isfile(cachefile), delete(cachefile); end
    warning('llm:i_cellcontext:BadDictionary', ...
        ['The cached marker dictionary was unreadable and has been ', ...
        'removed (%s). Annotation will run without the reference list.'], ...
        ME.message);
    recs = [];
    return
end

cached = recs;
end


function recs = i_validate(recs)
% Reject anything that is not the shape we expect. This file is fetched
% from the network and its contents are pasted into a model prompt, so it
% is treated as data of unverified shape, not as trusted input.
if ~isstruct(recs) || isempty(recs)
    error('llm:i_cellcontext:Shape', 'Expected a non-empty array of records.');
end
need = {'Tissue', 'CellType', 'Markers'};
if ~all(ismember(need, fieldnames(recs)))
    error('llm:i_cellcontext:Shape', ...
        'Records must have Tissue, Cell Type and Markers fields.');
end
ok = arrayfun(@(r) (ischar(r.Tissue) || isstring(r.Tissue)) && ...
    (ischar(r.CellType) || isstring(r.CellType)) && ...
    iscell(r.Markers), recs);
recs = recs(ok);
if isempty(recs)
    error('llm:i_cellcontext:Shape', 'No usable records.');
end
end


function s = i_clean(txt)
% One entry, reduced to a single line of plain text.
%
% Everything here ends up inside a prompt, so newlines and control
% characters are removed: without that, a crafted entry could close the
% reference section and forge instructions after it. Keeping every entry on
% one line means it can only ever read as one list item.
if isempty(txt)
    s = "";
    return
end
s = string(txt);
if ~isscalar(s), s = strjoin(s(:)', " "); end
s = regexprep(s, '[^\x20-\x7E]', ' ');          % printable ASCII only
s = regexprep(s, '[\[\]{}]', ' ');              % no bracket forging
s = strtrim(regexprep(s, '\s+', ' '));
s = extractBefore(s, min(strlength(s), 120) + 1);
end
