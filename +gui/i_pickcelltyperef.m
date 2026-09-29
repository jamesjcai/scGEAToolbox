function ref = i_pickcelltyperef(FigureHandle, source, querygenes)
%I_PICKCELLTYPEREF Ask for a cell-type reference and return it, or [] if none.
%
%   ref = gui.i_pickcelltyperef(FigureHandle)
%   ref = gui.i_pickcelltyperef(FigureHandle, "allen", sce.g)
%
% SOURCE "" (default) asks where the reference is; "allen" goes straight to
% the Allen Brain Cell Atlas. QUERYGENES, when given, lets the Allen route
% warn before matching a dataset that does not look like mouse against a
% mouse atlas: gene names are compared ignoring case, so a human dataset
% would otherwise run without complaint.
%
% The reference is the Allen Brain Cell Atlas (SC_ALLENBRAINREF, downloaded
% on first use after asking), or comes from a .mat file or the workspace,
% as either a struct made by SC_BUILDCELLTYPEREF (or assembled by hand from
% a published atlas's tables) or an annotated SCE, which is summarised on
% the spot and named after its file or variable. Every problem is reported
% in a dialog and returns [], so a caller only needs to check isempty.
%
% A list, not GUI.MYQUESTDLG: three sources plus the Cancel that
% MYQUESTDLG appends is one button more than it can show.
%
% see also: sc_allenbrainref, sc_buildcelltyperef,
% gui.callback_AnnotateByReference, gui.callback_AnnotationStability

arguments
    FigureHandle
    source (1, 1) string {mustBeMember(source, ["", "allen"])} = ""
    querygenes = []
end

ref = [];
if source == "allen"
    ref = in_allen(FigureHandle, querygenes);
    return;
end
sources = ["Allen Brain Cell Atlas (mouse brain)", ...
    "MAT file (reference or annotated SCE)", ...
    "Workspace variable (reference or annotated SCE)"];
[indx, tf] = gui.myListdlg(FigureHandle, sources, 'Select Reference', ...
    1, false, false, [340, 170], 'Where is the reference?');
if tf ~= 1, return; end
switch indx
    case 1
        ref = in_allen(FigureHandle, querygenes);
        return;
    case 2
        [value, sourcename] = in_fromfile(FigureHandle);
    case 3
        [value, sourcename] = in_fromworkspace(FigureHandle);
    otherwise
        % myListdlg with allowmulti=false returns one of the three indices.
        return;
end
if isempty(value), return; end

ref = in_toreference(FigureHandle, value, sourcename);
end


function ref = in_allen(FigureHandle, querygenes)
ref = [];
if ~isempty(querygenes) && ~strcmp(pkg.i_guessspecies(querygenes), 'mouse')
    msg = ['The gene names look human, and the Allen Brain Cell Atlas is ' ...
        'mouse. Genes are matched ignoring case, so this will run, but ' ...
        'it compares human cells with mouse cell types. Continue?'];
    if ~strcmp(gui.myQuestdlg(FigureHandle, msg, 'Allen Brain Cell Atlas', ...
            {'Continue', 'Cancel'}, 'Cancel'), 'Continue')
        return;
    end
end
try
    cached = sc_allenbrainref(AllowDownload=false);
catch ME
    if ~strcmp(ME.identifier, 'sc_allenbrainref:NotCached')
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    cached = [];
end
if isempty(cached)
    msg = ['The Allen Brain Cell Atlas reference is not on this computer ' ...
        'yet. Download it now from the Allen Institute (about 1.38 GB, ' ...
        'once; about 30 MB is kept)? The data is CC BY-NC 4.0: ' ...
        'non-commercial use only, with attribution to the Allen Institute ' ...
        '(Yao et al., Nature 2023).'];
    if ~strcmp(gui.myQuestdlg(FigureHandle, msg, 'Allen Brain Cell Atlas', ...
            {'Download', 'Cancel'}, 'Cancel'), 'Download')
        return;
    end
end

level = gui.myQuestdlg(FigureHandle, ['Match against the 338 subclasses ' ...
    'or the 34 classes of the whole mouse brain taxonomy?'], ...
    'Allen Brain Cell Atlas', {'Subclass', 'Class'}, 'Subclass');
if isempty(level) || ~ismember(level, {'Subclass', 'Class'}), return; end

fw = gui.myWaitbar(FigureHandle, [], false, ...
    'Preparing the Allen Brain Cell Atlas reference...');
try
    ref = sc_allenbrainref(Level=lower(string(level)));
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    ref = [];
    return;
end
gui.myWaitbar(FigureHandle, fw);
end


function [value, sourcename] = in_fromfile(FigureHandle)
value = [];
sourcename = "";
[fname, pathname] = uigetfile({'*.mat', 'MAT Files (*.mat)'}, ...
    'Select Reference or Annotated SCE File');
if isequal(fname, 0), return; end
scefile = fullfile(pathname, fname);

% Prefer a reference struct over an SCE: it is what a saved reference is,
% and building one from an SCE is the slower route.
info = whos('-file', scefile);
isRef = strcmp({info.class}, 'struct');
isSce = strcmp({info.class}, 'SingleCellExperiment');
if any(isRef)
    names = {info(isRef).name};
elseif any(isSce)
    names = {info(isSce).name};
else
    gui.myWarndlg(FigureHandle, sprintf(['%s holds neither a reference ' ...
        'struct nor a SingleCellExperiment.'], fname));
    return;
end
fw = gui.myWaitbar(FigureHandle);
try
    S = load(scefile, names{1});
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);
value = S.(names{1});
[~, sourcename] = fileparts(fname);
sourcename = string(sourcename);
end


function [value, sourcename] = in_fromworkspace(FigureHandle)
value = [];
sourcename = "";
a = evalin('base', 'whos');
isCandidate = ismember({a.class}, {'struct', 'SingleCellExperiment'});
if ~any(isCandidate)
    gui.myWarndlg(FigureHandle, ['No reference struct or SCE variable ' ...
        'in the workspace.']);
    return;
end
a = a(isCandidate);
labels = strcat({a.name}, ' (', {a.class}, ')');
[indx, tf] = gui.myListdlg(FigureHandle, labels, 'Select Reference', ...
    [], false, true, [], 'Pick the reference struct or annotated SCE.');
if tf ~= 1, return; end
value = evalin('base', a(indx).name);
sourcename = string(a(indx).name);
end


function ref = in_toreference(FigureHandle, value, sourcename)
ref = [];
if isa(value, 'SingleCellExperiment')
    if ~pkg.i_hascelltypelabels(value)
        gui.myWarndlg(FigureHandle, sprintf(['%s has no cell type ' ...
            'labels to use as a reference. Annotate it first.'], sourcename));
        return;
    end
    fw = gui.myWaitbar(FigureHandle);
    try
        ref = sc_buildcelltyperef(value.X, value.g, value.c_cell_type_tx, ...
            Source=sourcename);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        ref = [];
        return;
    end
    gui.myWaitbar(FigureHandle, fw);
elseif isstruct(value) && isscalar(value) && isfield(value, 'Types')
    ref = value;
    if ~isfield(ref, 'Source') || strlength(string(ref.Source)) == 0
        ref.Source = sourcename;
    end
else
    gui.myWarndlg(FigureHandle, sprintf(['%s is not a reference: it needs ' ...
        'the Types field of one made by sc_buildcelltyperef.'], sourcename));
end
end
