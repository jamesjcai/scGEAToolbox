function [ctype, colname] = i_pickcelltypecolumn(t, sce, parentfig, autoscore, ranked)
% I_PICKCELLTYPECOLUMN - Ask the user which metadata column holds cell types.
%
%   [ctype, colname] = gui.i_pickcelltypecolumn(t, sce)
%   [ctype, colname] = gui.i_pickcelltypecolumn(t, sce, parentfig, autoscore)
%   [ctype, colname] = gui.i_pickcelltypecolumn(t, sce, parentfig, autoscore, ranked)
%
% Inputs:
%   t         - per-cell metadata table, one row per cell
%   sce       - SingleCellExperiment the annotation will be applied to
%   parentfig - parent figure for the dialog (default: [])
%   autoscore - take the top ranked column without prompting when its score is
%               at least this value; Inf always prompts (default: Inf)
%   ranked    - ranking from a previous pkg.i_guesscelltypecol call on the same
%               table, to avoid scoring it twice; [] scores it here (default: [])
%
% Outputs:
%   ctype   - cell type string array with sce.NumCells elements, [] if cancelled
%   colname - name of the selected column, "" if cancelled
%
% Candidate columns are listed best-first using pkg.i_guesscelltypecol and are
% annotated with the number of distinct labels and a few example values, so an
% unfamiliar column can be judged from its contents rather than its name.
%
% see also: pkg.i_guesscelltypecol, gui.callback_AssignCellTypeFromAttrib

if nargin < 3, parentfig = []; end
if nargin < 4 || isempty(autoscore), autoscore = Inf; end
if nargin < 5, ranked = []; end

ctype = [];
colname = "";

if ~istable(t) || isempty(t), return; end
if ~isa(sce, 'SingleCellExperiment') || height(t) ~= sce.NumCells, return; end

if isempty(ranked)
    [bestname, bestscore, ranked] = pkg.i_guesscelltypecol(t, sce.NumCells);
else
    bestname = ranked.Name(1);
    bestscore = ranked.Score(1);
end
if isempty(ranked) || height(ranked) == 0, return; end

if bestscore >= autoscore
    colname = bestname;
else
    items = cell(height(ranked), 1);
    for k = 1:height(ranked)
        items{k} = sprintf('%s (%d types): %s', ranked.Name(k), ...
            ranked.NumTypes(k), ranked.Examples(k));
    end
    prompt = 'Select the column containing cell types:';
    if gui.i_isuifig(parentfig)
        [indx, tf] = gui.myListdlg(parentfig, items, prompt, items(1), false);
    else
        [indx, tf] = listdlg('PromptString', {prompt}, ...
            'SelectionMode', 'single', 'ListString', items, ...
            'InitialValue', 1, 'ListSize', [380, 300]);
    end
    if tf ~= 1 || isempty(indx), return; end
    colname = ranked.Name(indx);
end

ctype = strtrim(string(t.(colname)));
ctype(ismissing(ctype) | strlength(ctype) == 0) = "undetermined";
end
