function [k, genemode] = i_gethvgnum(sce, parentfig)
%I_GETHVGNUM Ask which genes an analysis should run on.
%
%   [k, genemode] = GUI.I_GETHVGNUM(sce, parentfig) returns the number of
%   genes K and a mode name GENEMODE, one of:
%
%     "hvg"                   the top K highly variable genes
%     "hvg+markers"           those plus every PanglaoDB marker present in
%                             the data
%     "hvg+markers:weighted"  the same set, with the markers the HVG cut
%                             had dropped scaled up so they carry more of
%                             the distance. Not supported by PHATE or
%                             METAVIZ, which normalise internally.
%     "all"                   every gene (K is then SCE.NUMGENES)
%
%   GENEMODE is "" and K is empty if the dialog was dismissed; callers must
%   test K and return.
%
%   GENEMODE passes straight to SINGLECELLEXPERIMENT.EMBEDCELLS's USEHVGS
%   argument, which reads these names and assembles the gene set itself.
%   Callers that select genes on their own instead take the marker list
%   from PKG.I_GETMARKERWHITELIST and append it with PKG.I_APPENDGENES.
%
%   More than three choices is too many for GUI.MYQUESTDLG: it reaches
%   QUESTDLG, which has no four-button syntax, or UICONFIRM, which caps at
%   four options and gets a fifth when MYQUESTDLG appends its own Cancel.
%   Hence a list dialog, which also has room to say what each choice costs.
%
%   See also PKG.I_GETMARKERWHITELIST, PKG.I_APPENDGENES.

if nargin < 2, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

k = [];
genemode = "";

defaulthvg = 2000;

% Count the markers the data actually contains, so the choice is not a leap
% of faith. This intersects the gene list only - no pass over the counts -
% so it stays cheap even on a 30000-gene SCE (~30 ms).
%
% The label reads as a union rather than a sum on purpose: most of these
% markers rank as highly variable anyway and are already in the top 2000,
% so the totals do not add. On the bundled 8260-cell dataset 840 of 1479
% overlap and the gene set comes to 2639, not 3479. Working out the true
% total needs the HVG ranking, which is the analysis itself, so the honest
% thing to show here is the two operands. EMBEDCELLS prints the number of
% genes it actually appended.
nmarker = numel(pkg.i_getmarkerwhitelist(sce.g));

items = { ...
    sprintf('%d HVGs', defaulthvg), ...
    sprintf('%d HVGs ∪ %d PanglaoDB markers 🐢', defaulthvg, nmarker), ...
    sprintf('%d HVGs ∪ %d PanglaoDB markers, markers up-weighted 🐢', ...
        defaulthvg, nmarker), ...
    sprintf('All genes (n=%d) 🐌', sce.NumGenes), ...
    'Other number of HVGs...'};

prompt = ['Which genes should the analysis use? Highly variable genes ' ...
    'are the fast default, but markers of rare cell types are often too ' ...
    'sparsely expressed to rank as variable and drop out. Adding the ' ...
    'PanglaoDB markers puts them back.'];

[indx, tf] = gui.myListdlg(parentfig, items, 'Select Genes', ...
    items{1}, false, true, [460, 190], prompt);
if ~tf, return; end

switch indx
    case 1
        genemode = "hvg";
        k = defaulthvg;
    case 2
        genemode = "hvg+markers";
        k = defaulthvg;
    case 3
        genemode = "hvg+markers:weighted";
        k = defaulthvg;
    case 4
        genemode = "all";
        k = sce.NumGenes;
    case 5
        k = gui.i_inputnumk(min([3000, sce.NumGenes]), ...
            100, sce.NumGenes, [], parentfig);
        if isempty(k), return; end
        genemode = "hvg";
    otherwise
        % MYLISTDLG returned TF true, so INDX is one of the five items.
        return;
end
end
