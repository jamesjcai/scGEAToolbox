function [Tm] = i_getcustomsubtypemarkers(parentfig, sce, targetname)
%I_GETCUSTOMSUBTYPEMARKERS Ask the user for a subtype marker list.
%
%   Tm = gui.i_getcustomsubtypemarkers(parentfig, sce, targetname)
%
%   Returns a two-column table - Var1 the subtype name, Var2 its markers as one
%   upper case comma separated string - or [] when the user cancels.
%
%   The list can be typed, loaded from a file, or started from the subtypes
%   assets/PanglaoDB/cellsubtypes.xlsx already has for TARGETNAME. Whichever it
%   is, it ends up in the editor before it is used, so a loaded file can still
%   be corrected and a typed list can still be saved by copying it out.
%
%   TARGETNAME is the label whose cells are about to be subdivided. It is used
%   only to pick a sensible starting point and to title the editor.
%
%   See also gui.callback_SubtypeAnnotation, pkg.i_parsemarkerlist,
%   pkg.i_readmarkertable, pkg.i_markerweights.

Tm = [];
if nargin < 3, targetname = ""; end
if nargin < 2, sce = []; end
if nargin < 1, parentfig = []; end

bundled = in_bundledsubtypes(targetname);

options = {'Type or paste a list...', 'Load from a file...'};
if strlength(bundled) > 0
    options = [{'Start from the bundled subtypes...'}, options];
end
answer = gui.myQuestdlg(parentfig, ...
    sprintf(['Subtype markers for %s.\n\nFormat: one subtype per line, ' ...
    'name and genes separated by a tab.'], in_display(targetname)), ...
    'Customized Subtype Markers', [options, {'Cancel'}], options{1});
if isempty(answer) || strcmp(answer, 'Cancel'), return; end

switch answer
    case 'Start from the bundled subtypes...'
        indata = bundled;
    case 'Load from a file...'
        [f, p] = uigetfile({'*.txt;*.csv;*.tsv;*.text', ...
            'Marker list (*.txt, *.csv, *.tsv)'; ...
            '*.xlsx;*.xls', 'Excel marker table (*.xlsx, *.xls)'; ...
            '*.*', 'All files (*.*)'}, 'Select a marker gene file');
        if isequal(f, 0), return; end
        try
            indata = pkg.i_readmarkertable(fullfile(p, f));
        catch ME
            gui.myErrordlg(parentfig, ME.message, 'Marker file');
            return;
        end
        if strlength(indata) == 0
            gui.myErrordlg(parentfig, sprintf(['No marker list could be ' ...
                'read from %s. Expected one cell type per line, its name ' ...
                'and its genes separated by a tab, or a spreadsheet whose ' ...
                'first two columns hold the same.'], f), 'Marker file');
            return;
        end
    otherwise
        % Typed from scratch: a template built from genes that are actually in
        % this dataset, so the format is unambiguous and the example is real.
        indata = in_template(sce, targetname);
end

if gui.i_isuifig(parentfig)
    a = gui.myInputwin([], [], char(indata), parentfig);
else
    a = inputdlg(sprintf('Format:\nSubtype name [TAB] Gene1,Gene2'), ...
        'Subtype Markers Input', [15, 80], {char(indata)}, 'on');
end
if isempty(a), return; end

Tm = pkg.i_parsemarkerlist(a);
if isempty(Tm)
    gui.myErrordlg(parentfig, ['No usable marker list. Each line needs a ' ...
        'subtype name, a tab, and at least one gene symbol.'], ...
        'Customized Subtype Markers');
    Tm = [];
    return;
end
if height(Tm) < 2
    % One subtype is not a partition: every cluster scores against the same
    % list and the whole population comes back with one label.
    answer = gui.myQuestdlg(parentfig, sprintf(['Only one subtype (%s) was ' ...
        'given, so every cell of %s will end up with that label. ' ...
        'Continue?'], Tm.Var1(1), in_display(targetname)), ...
        'Customized Subtype Markers', {'Continue', 'Cancel'}, 'Cancel');
    if ~strcmp(answer, 'Continue'), Tm = []; return; end
end
end

function [txt] = in_bundledsubtypes(targetname)
% The subtypes cellsubtypes.xlsx lists for the primary type TARGETNAME belongs
% to, as editor text. "" when the label reaches no primary type, which is the
% usual case for the cell types this menu item exists for.

txt = "";
if strlength(targetname) == 0, return; end

pw1 = fileparts(fileparts(mfilename('fullpath')));
pth = fullfile(pw1, 'assets', 'PanglaoDB', 'cellsubtypes.xlsx');
if ~exist(pth, 'file'), pth = which('cellsubtypes.xlsx'); end
if isempty(pth) || ~exist(pth, 'file'), return; end

try
    T = readtable(pth, 'TextType', 'string');
catch
    return;
end
if ~all(ismember({'CellType', 'SubType', 'PositiveMarkers'}, ...
        T.Properties.VariableNames))
    return;
end

primary = pkg.i_matchprimarytype(targetname, unique(T.CellType));
if strlength(primary) == 0, return; end

T = T(upper(T.CellType) == upper(primary), :);
if isempty(T), return; end
txt = strjoin(strtrim(T.SubType) + sprintf('\t') + ...
    upper(erase(strtrim(T.PositiveMarkers), " ")), newline);
end

function [txt] = in_template(sce, targetname)
% A two-line example. Real gene symbols from the data beat "gene1,gene2":
% they show the casing the dataset uses, which is what decides whether a
% hand-typed marker matches anything at all.

name = in_display(targetname);
if isempty(sce) || ~isa(sce, 'SingleCellExperiment') || numel(sce.g) < 7
    g = "GENE" + (1:7)';
else
    rngState = rng();
    restoreRng = onCleanup(@() rng(rngState));
    rng("shuffle");
    g = string(sce.g(randperm(numel(sce.g), 7)));
end

txt = sprintf('%s subtype 1\t%s,%s,%s,%s\n%s subtype 2\t%s,%s,%s', ...
    name, g(1), g(2), g(3), g(4), name, g(5), g(6), g(7));
end

function [s] = in_display(targetname)
s = strtrim(string(targetname));
if strlength(s) == 0, s = "the selected cells"; end
end
