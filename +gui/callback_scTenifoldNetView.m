function callback_scTenifoldNetView(src, ~)
% CALLBACK_SCTENIFOLDNETVIEW  View a saved scTenifoldNet result.
%
% Opens the side-by-side viewer for a pair of GRNs produced by
% ten.sctenifoldnet or gui.callback_scTenifoldNet1lite, with an optional
% third panel for their difference. The networks come either from two .mat
% files on disk or from two matrices already in the base workspace; the
% differential regulation table, if there is one, decides which genes are
% drawn.
%
% See also ten.sctenifoldnetview, gui.callback_scTenifoldNet2.

if nargin < 1, src = []; end
[FigureHandle, ~] = gui.gui_getfigsce(src);

if ~gui.gui_showrefinfo('scTenifoldNet [PMID:33336197]', FigureHandle), return; end

answer = gui.myQuestdlg(FigureHandle, ...
    'Where are the two networks?', 'Input Networks', ...
    {'Files (.mat)', 'Workspace', 'Cancel'}, 'Files (.mat)');

switch answer
    case 'Files (.mat)'
        [A0, A1, genes, T, labels, celltype] = i_fromfiles(FigureHandle);
    case 'Workspace'
        [A0, A1, genes, T, labels, celltype] = i_fromworkspace(FigureHandle);
    otherwise
        return;
end
if isempty(A0) || isempty(A1), return; end

modelist = {'Top DR genes and their partners', ...
    'All significant DR genes (the whole DR module)'};
[modeidx, ok] = gui.myListdlg(FigureHandle, modelist, 'Which Genes', ...
    [], false, true, [340, 200], ...
    ['Pick the gene set to draw. The second option keeps every gene ', ...
    'below the significance cutoff instead of a top-N slice, so the ', ...
    'number of seeds below is ignored.']);
if ~ok, return; end
if modeidx == 2
    viewmode = "drmodule";
else
    viewmode = "neighborhood";
end

para = gui.myInputdlg({'Number of seed (top DR) genes:', ...
    'Total number of nodes to draw:', ...
    'Edge cutoff quantile (0-1, higher = fewer edges):', ...
    'Minimum edges kept per node:'}, ...
    'Viewer Settings', {'8', '60', '0.9', '2'}, FigureHandle);
if isempty(para), return; end
numseeds = str2double(para{1});
numnodes = str2double(para{2});
edgecutoff = str2double(para{3});
mindegree = str2double(para{4});
if ~all(isfinite([numseeds, numnodes, edgecutoff, mindegree])) ...
        || numseeds < 1 || numnodes < 2 || edgecutoff < 0 || edgecutoff >= 1 ...
        || mindegree < 0
    gui.myErrordlg(FigureHandle, ['Settings must be numeric: at least 1 seed, ', ...
        'at least 2 nodes, an edge cutoff in [0, 1), and a non-negative ', ...
        'minimum degree.']);
    return;
end

showdiff = strcmp('Yes', gui.myQuestdlg(FigureHandle, ...
    ['Add a third panel showing the difference between the two networks? ', ...
    'It shares their node layout, so an edge that changed is read off ', ...
    'the same position.'], 'Difference Panel', {'Yes', 'No'}, 'No'));

fw = gui.myWaitbar(FigureHandle);
try
    ten.sctenifoldnetview(A0, A1, genes, T, ...
        'Mode', viewmode, ...
        'NumSeeds', round(numseeds), ...
        'NumNodes', round(numnodes), ...
        'EdgeCutoff', edgecutoff, ...
        'MinDegree', round(mindegree), ...
        'ShowDifference', showdiff, ...
        'Labels', labels, ...
        'Title', celltype, ...
        'ParentFig', FigureHandle);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
end

% ------------------------------------------------------- local functions

function [A0, A1, genes, T, labels, celltype] = i_fromfiles(parentfig)
[A0, A1, genes, T, labels, celltype] = deal([], [], [], [], string.empty, string.empty);

[f1, p1] = uigetfile({'*.mat'}, 'Select the FIRST network (.mat)');
if isequal(f1, 0), return; end
[f2, p2] = uigetfile({'*.mat'}, 'Select the SECOND network (.mat)', p1);
if isequal(f2, 0), return; end

% ten.sctenifoldnetview reads A, the gene list and the metadata itself.
A0 = fullfile(p1, f1);
A1 = fullfile(p2, f2);

[labels, celltype] = i_readlabels(A0, A1);

if strcmp('Yes', gui.myQuestdlg(parentfig, ...
        ['Load a differential regulation (DR) table? Without one, genes ', ...
        'are ranked by how far their edges moved instead.'], ...
        'DR Table'))
    [f3, p3] = uigetfile({'*.csv;*.txt', 'Table (*.csv, *.txt)'; '*.mat', 'MAT-file (*.mat)'}, ...
        'Select the DR table', p1);
    if ~isequal(f3, 0)
        T = i_readdrtable(fullfile(p3, f3), parentfig);
    end
end
end

function [labels, celltype] = i_readlabels(f1, f2)
labels = string.empty;
celltype = string.empty;
try
    m1 = load(f1, 'genotype', 'celltype');
    m2 = load(f2, 'genotype');
    if isfield(m1, 'genotype') && isfield(m2, 'genotype')
        labels = [string(m1.genotype), string(m2.genotype)];
    end
    if isfield(m1, 'celltype')
        celltype = string(m1.celltype);
    end
catch
    % Metadata is optional; the viewer falls back to generic labels.
end
end

function T = i_readdrtable(fname, parentfig)
T = [];
try
    [~, ~, ext] = fileparts(fname);
    if strcmpi(ext, '.mat')
        S = load(fname);
        fn = fieldnames(S);
        istbl = cellfun(@(x) istable(S.(x)), fn);
        if ~any(istbl)
            gui.myErrordlg(parentfig, 'That MAT-file holds no table.');
            return;
        end
        T = S.(fn{find(istbl, 1)});
    else
        T = readtable(fname);
    end
catch ME
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    T = [];
end
end

function [A0, A1, genes, T, labels, celltype] = i_fromworkspace(parentfig)
[A0, A1, genes, T, labels, celltype] = deal([], [], [], [], string.empty, string.empty);

if ismcc || isdeployed
    gui.myErrordlg(parentfig, 'Workspace access is not available in the standalone application.');
    return;
end
vars = evalin('base', 'whos');
if isempty(vars)
    gui.myHelpdlg(parentfig, 'No variable in the Workspace.');
    return;
end
varnames = {vars.name};

% Two single-selection dialogs rather than one multi-select: which network
% is the reference matters, and a multi-select returns list order, not the
% order the user clicked in.
[i1, tf] = gui.myListdlg(parentfig, varnames, 'First Network', [], false, true, ...
    [300, 450], 'Select the adjacency matrix of the FIRST (reference) condition.');
if ~tf, return; end
[i2, tf] = gui.myListdlg(parentfig, varnames, 'Second Network', [], false, true, ...
    [300, 450], 'Select the adjacency matrix of the SECOND condition.');
if ~tf, return; end
A0 = evalin('base', varnames{i1});
A1 = evalin('base', varnames{i2});
labels = string(varnames([i1, i2]));
if ~isequal(size(A0), size(A1))
    gui.myErrordlg(parentfig, 'The two networks must be the same size.');
    [A0, A1] = deal([]);
    return;
end

[indx, tf] = gui.myListdlg(parentfig, varnames, 'Gene List', [], false, true, ...
    [300, 450], 'Select the gene list indexing both networks.');
if ~tf
    [A0, A1] = deal([]);
    return;
end
genes = string(evalin('base', varnames{indx}));
if numel(genes) ~= size(A0, 1)
    gui.myErrordlg(parentfig, sprintf( ...
        'The gene list has %d entries but the networks are %dx%d.', ...
        numel(genes), size(A0, 1), size(A0, 2)));
    [A0, A1] = deal([]);
    return;
end

if strcmp('Yes', gui.myQuestdlg(parentfig, ...
        ['Use a differential regulation (DR) table from the Workspace? ', ...
        'Without one, genes are ranked by how far their edges moved instead.'], ...
        'DR Table'))
    istbl = arrayfun(@(v) strcmp(v.class, 'table'), vars);
    if ~any(istbl)
        gui.myHelpdlg(parentfig, 'No table in the Workspace.');
        return;
    end
    tblnames = varnames(istbl);
    [indx, tf] = gui.myListdlg(parentfig, tblnames, 'DR Table', [], false, true, ...
        [300, 450], 'Select the DR table returned by ten.sctenifoldnet.');
    if tf
        T = evalin('base', tblnames{indx});
    end
end
end
