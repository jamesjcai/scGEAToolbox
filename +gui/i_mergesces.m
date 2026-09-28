function sce = i_mergesces(parentfig, insce, names, confirmfcn)
% I_MERGESCES - Ask how to merge several SCEs, then merge them.
%
%   sce = gui.i_mergesces(parentfig, insce, names)
%   sce = gui.i_mergesces(parentfig, insce, names, confirmfcn)
%
% INSCE is a cell array of two or more SingleCellExperiment objects and NAMES
% a string array naming each one (a variable or file name). NAMES label the
% cells' batches when the inputs' own batch IDs cannot tell them apart.
%
% Asks at most two questions, each only when it has more than one sensible
% answer: which genes to keep, when the gene lists differ; and how to label
% batches, when some input already holds several. CONFIRMFCN, when given, is
% called with the number of merged cells once both are answered and before
% any work is done; returning false cancels.
%
% Returns [] when cancelled or when the merge fails, having said why.
%
% This is the one merge dialog behind File > Import Data (several SCE files,
% or several workspace variables) and Edit > Merge Current Dataset with
% Others.
%
% see also: sc_mergesces, gui.i_loadscefiles, gui.i_pickworkspacesces

if nargin < 4, confirmfcn = []; end
sce = [];
names = matlab.lang.makeUniqueStrings(string(names(:)'));

methodtag = in_pickgenemethod(parentfig, insce);
if isempty(methodtag), return; end

batchmode = in_pickbatchmode(parentfig, insce, names);
if isempty(batchmode), return; end

if ~isempty(confirmfcn)
    ncells = sum(cellfun(@(x) x.NumCells, insce));
    if ~confirmfcn(ncells), return; end
end

fw = gui.myWaitbar(parentfig);
try
    % Batch labels are set here rather than left to SC_MERGESCES, which
    % numbers the inputs 1..n and appends those numbers on a collision; the
    % names say which input a cell came from.
    merged = sc_mergesces(insce, methodtag, true, true);
    merged.c_batch_id = in_batchlabels(insce, names, batchmode);
catch ME
    gui.myWaitbar(parentfig, fw, true);
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(parentfig, fw);
sce = merged;
end

function methodtag = in_pickgenemethod(parentfig, insce)
% '' when cancelled, or when the inputs share no gene at all.

gi = string(insce{1}.g(:));
gu = gi;
for k = 2:numel(insce)
    gk = string(insce{k}.g(:));
    gi = intersect(gi, gk);
    gu = union(gu, gk);
end

methodtag = '';
if isempty(gi)
    gui.myWarndlg(parentfig, ['The selected datasets have no gene names in ' ...
        'common, so they cannot be merged. They may come from different ' ...
        'species, or use different gene identifiers (symbols in one, ' ...
        'Ensembl IDs in another).']);
    return;
end
if numel(gi) == numel(gu)
    methodtag = 'intersect';
    return;
end

msg = sprintf(['The datasets share %d genes; together they have %d.\n\n' ...
    'Shared genes only: every gene kept is measured in every cell.\n' ...
    'All genes: a gene missing from one dataset is set to zero in its ' ...
    'cells.'], numel(gi), numel(gu));
answer = gui.myQuestdlg(parentfig, msg, 'Genes to Keep', ...
    {'Shared Genes Only', 'All Genes'}, 'Shared Genes Only');
switch answer
    case 'Shared Genes Only'
        methodtag = 'intersect';
    case 'All Genes'
        methodtag = 'union';
    otherwise
        % Dismissed; methodtag stays ''.
end
end

function batchmode = in_pickbatchmode(parentfig, insce, names)
% 'original' keeps each cell's batch ID, 'source' labels every cell with the
% input it came from, '' means cancelled.
%
% Only asked when some input already holds several batches. Otherwise each
% input has one batch ID at most: those are kept when they tell the inputs
% apart (sample accessions, say), and replaced by the names when they do not
% (the default "1" in every input).

n = numel(insce);
batches = cell(1, n);
for k = 1:n
    batches{k} = unique(in_batchof(insce{k}));
end
nb = cellfun(@numel, batches);

if all(nb <= 1)
    orig = vertcat(batches{:});
    if all(nb == 1) && numel(unique(orig)) == n
        batchmode = 'original';
    else
        batchmode = 'source';
    end
    return;
end

multi = names(nb > 1);
msg = sprintf(['%s already contains several batches. Label the merged ' ...
    'cells by their original batch IDs (%d batches), or by the dataset ' ...
    'each came from (%d batches: %s)?\n\nWhere two datasets use the same ' ...
    'batch ID, the original IDs are prefixed with the dataset name.'], ...
    strjoin(multi, ", "), sum(nb), n, strjoin(names, ", "));
answer = gui.myQuestdlg(parentfig, msg, 'Batch Labels', ...
    {'Original Batch IDs', 'One per Dataset'}, 'Original Batch IDs');
switch answer
    case 'Original Batch IDs'
        batchmode = 'original';
    case 'One per Dataset'
        batchmode = 'source';
    otherwise
        batchmode = '';
end
end

function b = in_batchof(sce)
% Batch IDs as a string column, or empty when none of the right length.
b = string(sce.c_batch_id(:));
if numel(b) ~= sce.NumCells
    b = strings(0, 1);
end
end

function b = in_batchlabels(insce, names, batchmode)
% One label per merged cell, in the order SC_MERGESCES concatenates them.

n = numel(insce);
parts = cell(n, 1);
for k = 1:n
    nk = insce{k}.NumCells;
    bk = in_batchof(insce{k});
    if strcmp(batchmode, 'source') || isempty(bk)
        parts{k} = repmat(names(k), nk, 1);
    else
        parts{k} = bk;
    end
end

if strcmp(batchmode, 'original')
    % A batch ID used by two inputs would put unrelated cells in one batch,
    % so in that case every ID is prefixed with its input's name.
    allids = vertcat(parts{:});
    owner = repelem((1:n)', cellfun(@numel, parts));
    [~, ~, g] = unique(allids);
    nowners = accumarray(g, owner, [], @(x) numel(unique(x)));
    if any(nowners > 1)
        for k = 1:n
            parts{k} = names(k) + "_" + parts{k};
        end
    end
end
b = vertcat(parts{:});
end
