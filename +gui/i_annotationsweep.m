function [S, args, speciestag] = i_annotationsweep(src, FigureHandle, sce)
%I_ANNOTATIONSWEEP The resolution sweep of SC_ANNOTATIONSTABILITY, run once per dataset.
%
%   [S, args, speciestag] = gui.i_annotationsweep(app, FigureHandle, sce)
%
% The sweep clusters and labels the data at every resolution, which is the
% slow part of both menu items that use it (Cluster Cells > Louvain on PCs >
% Choose by Annotation, and Cluster > Choose Clustering Resolution by
% Annotation). Its result is kept in sce.struct_saved_results with a
% fingerprint of the data, and while the data still match it is offered
% again instead of being recomputed. Changing the cells, genes or counts, or
% running Harmony, changes the fingerprint and the saved result is dropped.
%
% Either way the per-resolution table is shown, with a Methods Text button
% for the paragraph SC_ANNOTATIONSTABILITY writes.
%
% S is [] when the user cancels, otherwise a struct:
%   Method, MethodLabel   the labelling method, and how to name it
%   Species, Markers      the options it used ("markers", "custommarkers")
%   ReferenceSource       what the reference was ("reference")
%   Date                  when it ran
%   T                     SC_ANNOTATIONSTABILITY's table; its
%                         Properties.Description is the Methods paragraph
%   Resolution, ClusterId the recommended resolution and its partition
%   Stability             each cell's share of resolutions agreeing on its type
%   Fingerprint           of the data it was computed on
%
% ARGS are the name-value arguments that label clusters the same way through
% SC_ANNOTATECELLS. They are {} for a saved reference sweep: the reference
% is not kept with the data (the Allen atlas alone is tens of megabytes), so
% a caller that needs to label must ask for it again. SPECIESTAG is the
% species chosen for "markers", [] otherwise.
%
% see also: sc_annotationstability, gui.callback_AnnotationStability,
% gui.callback_ReclusterCells

S = [];
args = {};
speciestag = [];

fingerprint = i_fingerprint(sce);
saved = [];
if isfield(sce.struct_saved_results, 'annotationstability')
    saved = sce.struct_saved_results.annotationstability;
    if ~isequal(saved.Fingerprint, fingerprint)
        % The data changed since the sweep ran; it no longer describes them.
        sce.struct_saved_results = rmfield(sce.struct_saved_results, ...
            'annotationstability');
        saved = [];
    end
end

if ~isempty(saved)
    msg = sprintf(['The resolution sweep was already run on these data ' ...
        '(%s, %s) and recommended resolution %g. Use that result, or run ' ...
        'the sweep again?'], saved.MethodLabel, ...
        string(saved.Date, 'yyyy-MM-dd HH:mm'), saved.Resolution);
    switch gui.myQuestdlg(FigureHandle, msg, 'Choose Clustering Resolution', ...
            {'Use Saved Result', 'Run Again'}, 'Use Saved Result')
        case 'Use Saved Result'
            S = saved;
            i_showtable(FigureHandle, S.T);
            switch S.Method
                case "markers"
                    speciestag = char(S.Species);
                    args = {'Species', S.Species};
                case "custommarkers"
                    args = {'Markers', S.Markers};
                otherwise
                    % "reference": not kept with the data; see help.
                    args = {};
            end
            return;
        case 'Run Again'
            % Fall through to the method choice below.
        otherwise
            % Cancel, or the dialog closed.
            return;
    end
end

% A list, not GUI.MYQUESTDLG: three methods plus the Cancel it appends is
% one button more than it can show.
methods = ["Database markers (PanglaoDB)", "Custom marker list", ...
    "Reference (Allen Brain Cell Atlas or your own)"];
[indx, tf] = gui.myListdlg(FigureHandle, methods, ...
    'Choose Clustering Resolution', 1, false, false, [360, 170], ...
    ['Each resolution''s clusters are labelled, and the coarsest ' ...
    'resolution that shows every cell type found reliably is ' ...
    'recommended. Label them with:']);
if tf ~= 1, return; end
species = "";
markers = table();
refSource = "";
switch indx
    case 1
        method = "markers";
        preferred = [];
        if isprop(src, 'speciestag'), preferred = src.speciestag; end
        speciestag = gui.i_selectspecies(2, false, FigureHandle, preferred);
        if isempty(speciestag), return; end
        species = string(speciestag);
        args = {'Species', species};
        label = sprintf('PanglaoDB markers, %s', species);
    case 2
        method = "custommarkers";
        markers = gui.i_getcustommarkers(FigureHandle, sce);
        if isempty(markers), return; end
        args = {'Markers', markers};
        label = sprintf('custom list of %d cell types', height(markers));
    case 3
        method = "reference";
        ref = gui.i_pickcelltyperef(FigureHandle, "", sce.g);
        if isempty(ref), return; end
        args = {'Reference', ref};
        if isfield(ref, 'Source') && strlength(string(ref.Source)) > 0
            refSource = string(ref.Source);
        end
        label = sprintf('reference of %d cell types', numel(ref.Types));
        if refSource ~= ""
            label = sprintf('%s (%s)', label, refSource);
        end
    otherwise
        % myListdlg with allowmulti=false returns one of the three indices.
        return;
end

fw = gui.myWaitbar(FigureHandle, [], false, ...
    'Clustering and labelling at each resolution...');
try
    [T, bestResolution, bestClusterId, perCell] = sc_annotationstability( ...
        sce, 'Method', method, args{:}, 'Verbose', false);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    args = {};
    speciestag = [];
    return;
end
gui.myWaitbar(FigureHandle, fw);

S = struct('Method', method, 'MethodLabel', string(label), ...
    'Species', species, 'Markers', markers, 'ReferenceSource', refSource, ...
    'Date', datetime('now'), 'T', T, 'Resolution', bestResolution, ...
    'ClusterId', bestClusterId, 'Stability', perCell.Stability, ...
    'Fingerprint', fingerprint);
sce.struct_saved_results.annotationstability = S;
i_showtable(FigureHandle, T);
end


function i_showtable(FigureHandle, T)
% The per-resolution table, with the Methods paragraph a button away.
methodsAction = struct('Text', 'Methods Text', ...
    'Tooltip', 'A paragraph describing this run, for the Methods section of a paper', ...
    'Callback', @(~, viewerFig) gui.myTextareadlg(viewerFig, {'Methods:'}, ...
    'Methods Text', {char(T.Properties.Description)}, true));
gui.TableViewerApp(T, FigureHandle, "ResolutionSweep", methodsAction);
end


function f = i_fingerprint(sce)
% Cheap, order-sensitive summary of what the sweep depends on: the counts,
% the gene names, and whether Harmony components replace the PCs.
X = sce.X;
[numGenes, numCells] = size(X);
cellSums = full(sum(X, 1));
geneSums = full(sum(X, 2));
r = sce.struct_cell_reductions;
harmony = [];
if isstruct(r) && isfield(r, 'harmony') && ~isempty(r.harmony)
    harmony = [size(r.harmony), double(sum(r.harmony, 'all'))];
end
f = struct('Size', [numGenes, numCells], 'Nnz', nnz(X), ...
    'CellSums', double(cellSums*(1:numCells)'), ...
    'GeneSums', double((1:numGenes)*geneSums), ...
    'Total', double(sum(cellSums)), ...
    'Genes', string(pkg.i_toHash(char(strjoin(sce.g, newline)))), ...
    'Harmony', harmony);
end
