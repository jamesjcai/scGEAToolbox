function [Tm, srcname] = i_getcustommarkers(parentfig, sce, dlgtitle)
%I_GETCUSTOMMARKERS Ask the user for a cell type marker list.
%
%   [Tm, srcname] = gui.i_getcustommarkers(parentfig, sce)
%   [Tm, srcname] = gui.i_getcustommarkers(parentfig, sce, dlgtitle)
%
%   Returns a two-column table - Var1 the cell type name, Var2 its markers as
%   one upper case comma separated string - or [] when the user cancels.
%   SRCNAME says where the list came from, for a caller that wants to name it
%   in a later message.
%
%   The list can be typed, loaded from a file, or started from the bundled
%   ScTypeDB markers for a tissue. Whichever it is, it ends up in the editor
%   before it is used, so a loaded list can still be corrected.
%
%   Every usable list is remembered for the rest of the MATLAB session through
%   GUI.I_SESSIONMARKERS, and from the second run on it leads the menu and is
%   the selected row: pressing OK confirms the markers just used, and any
%   other row changes them. Typing a marker list is slow enough that retyping
%   it for the next selection is the main reason this path goes unused, and a
%   list that is only offered costs nothing to remember.
%
%   See also gui.i_sessionmarkers, gui.i_getcustomsubtypemarkers,
%   pkg.i_parsemarkerlist, pkg.i_readmarkertable, gui.i_warnmissingmarkers.

Tm = [];
srcname = "";
if nargin < 3 || isempty(dlgtitle), dlgtitle = 'Customized Marker Genes'; end
if nargin < 2, sce = []; end
if nargin < 1, parentfig = []; end

[savedTm, savedname] = gui.i_sessionmarkers();
hassaved = ~isempty(savedTm);

items = ["Type or paste a list..."; ...
    "Load from a file..."; ...
    "Start from the bundled ScTypeDB markers..."];
prompt = ['Where should the marker genes come from? A list is one cell ' ...
    'type per line, its name and its genes separated by a tab.'];

if hassaved
    % First, and so the selected row - a listbox opens on its first item -
    % which makes the common case, the same markers on the next selection,
    % one OK, and changing them one click rather than a retype.
    %
    % The prompt names them. "The saved list" is not something the user can
    % check against what they meant to use, and this is the last chance to
    % notice that the markers are from the run before last, or for the other
    % species, or the ones typed to test a hunch.
    items = ["Use these markers again"; "Edit these markers..."; items];
    prompt = sprintf(['The markers used last time are %d cell type(s) ' ...
        'from %s: %s. Press OK to use them again, or pick another row ' ...
        'to change them.'], height(savedTm), savedname, ...
        pkg.i_namesummary(savedTm.Var1, 3));
end

[indx, tf] = gui.myListdlg(parentfig, cellstr(items), dlgtitle, ...
    [], false, true, [420, 260], prompt);
if tf ~= 1 || isempty(indx), return; end
choice = string(items(indx));

if choice == "Use these markers again"
    % Already parsed, and already checked against the data when it arrived.
    Tm = savedTm;
    srcname = savedname;
    return;
end

switch choice
    case "Edit these markers..."
        indata = in_tabletext(savedTm);
        % Not the file's contents any more once the editor has been through
        % them, and the next run names the source back to the user.
        srcname = savedname;
        if ~endsWith(srcname, " (edited)")
            srcname = srcname + " (edited)";
        end
    case "Load from a file..."
        [f, p] = uigetfile({'*.txt;*.csv;*.tsv;*.text', ...
            'Marker list (*.txt, *.csv, *.tsv)'; ...
            '*.xlsx;*.xls', 'Excel marker table (*.xlsx, *.xls)'; ...
            '*.*', 'All files (*.*)'}, 'Select a marker gene file');
        if isequal(f, 0), return; end
        try
            indata = pkg.i_readmarkertable(fullfile(p, f));
        catch ME
            gui.myErrordlg(parentfig, ME.message, dlgtitle);
            return;
        end
        if strlength(indata) == 0
            gui.myErrordlg(parentfig, sprintf(['No marker list could be ' ...
                'read from %s. Expected one cell type per line, its name ' ...
                'and its genes separated by a tab, or a spreadsheet whose ' ...
                'first two columns hold the same.'], f), dlgtitle);
            return;
        end
        srcname = string(f);
    case "Start from the bundled ScTypeDB markers..."
        indata = string(gui.i_getsctypemarkers(parentfig));
        if strlength(indata) == 0, return; end
        srcname = "ScTypeDB";
    otherwise
        % Typed from scratch: a template built from genes that are actually in
        % this dataset, so the format is unambiguous and the example is real.
        indata = in_template(sce);
        srcname = "typed list";
end

if gui.i_isuifig(parentfig)
    a = gui.myInputwin([], [], char(indata), parentfig);
else
    a = inputdlg(sprintf('Format:\nCell type name [TAB] Gene1,Gene2'), ...
        dlgtitle, [15, 80], {char(indata)}, 'on');
end
if isempty(a), return; end

Tm = pkg.i_parsemarkerlist(a);
if isempty(Tm)
    gui.myErrordlg(parentfig, ['No usable marker list. Each line needs a ' ...
        'cell type name, a tab, and at least one gene symbol.'], dlgtitle);
    Tm = [];
    return;
end

gui.i_warnmissingmarkers(parentfig, sce, Tm, dlgtitle);

% Remembered without asking. The next run leads with it and says what it is,
% so the user confirms the markers where they are about to be used rather
% than answering a question about a list they have not seen since.
gui.i_sessionmarkers(Tm, srcname);
end

function [txt] = in_tabletext(Tm)
% A parsed marker table back as editor text, so a saved list can be corrected
% and extended in the same editor it was typed in.

txt = strjoin(strtrim(string(Tm.Var1)) + sprintf('\t') + ...
    strtrim(string(Tm.Var2)), newline);
end

function [txt] = in_template(sce)
% A two-line example. Real gene symbols from the data beat "gene1,gene2":
% they show the casing the dataset uses, which is what decides whether a
% hand-typed marker matches anything at all.

if isempty(sce) || ~isa(sce, 'SingleCellExperiment') || numel(sce.g) < 7
    g = "GENE" + (1:7)';
else
    rngState = rng();
    restoreRng = onCleanup(@() rng(rngState));
    rng("shuffle");
    g = string(sce.g(randperm(numel(sce.g), 7)));
end

txt = sprintf('Cell type 1\t%s,%s,%s,%s\nCell type 2\t%s,%s,%s', ...
    g(1), g(2), g(3), g(4), g(5), g(6), g(7));
end
