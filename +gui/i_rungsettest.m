function T = i_rungsettest(genes, stats, parentfig, defname, statlabel)
%I_RUNGSETTEST  Competitive gene-set test on a ranked gene list, with dialogs.
%
%   T = GUI.I_RUNGSETTEST(genes, stats, parentfig, defname) prompts for a
%   gene set collection and a test method, runs SC_GSETTEST on the ranking
%   defined by STATS, and shows the result in GUI.TABLEVIEWERAPP. T is empty
%   if the user cancels either prompt.
%
%   STATS must be the statistic for EVERY measured gene, not a filtered or
%   truncated top-N list. A competitive test asks whether a set's genes rank
%   higher than the genes outside it, so the genes that failed a significance
%   cut are exactly the comparison group. Handing this function a pre-filtered
%   marker list silently shrinks the universe to those markers, and the test
%   stops meaning anything.
%
%   Sets larger than 500 genes are skipped, following the usual GSEA/fgsea
%   convention: very large GO terms otherwise crowd the top of the result and
%   say little about biology.
%
%   INPUTS:
%     genes     - G-by-1 gene symbols
%     stats     - G-by-1 ranking statistic, larger meaning "more up"
%     parentfig - figure to centre the dialogs on
%     defname   - default file name offered by the table viewer's export
%     statlabel - name of the statistic, shown in the method prompt so the
%                 user can see what is being ranked (optional)
%
%   USAGE (as a GUI.TABLEVIEWERAPP custom action):
%     action = struct('Text', 'Gene-Set Test', 'Tooltip', '...', ...
%         'Callback', @(~, figtab) gui.i_rungsettest(string(T.gene), ...
%             T.avg_log2FC, figtab, 'DE_GSet', 'avg_log2FC'));
%
% See also: SC_GSETTEST, PKG.E_GETGENESETS, GUI.TABLEVIEWERAPP

if nargin < 5, statlabel = ''; end
if nargin < 4 || isempty(defname), defname = 'GeneSetTest'; end
if nargin < 3, parentfig = []; end

T = table.empty;

genes = string(genes(:));
stats = double(stats(:));
if numel(genes) ~= numel(stats)
    gui.myErrordlg(parentfig, sprintf( ...
        'Gene list (%d) and statistic (%d) must have the same length.', ...
        numel(genes), numel(stats)));
    return;
end

ok = isfinite(stats) & strlength(genes) > 0;
if sum(ok) < 10
    gui.myErrordlg(parentfig, sprintf( ...
        ['Only %d genes have a usable statistic. A competitive gene-set ' ...
        'test needs the whole measured gene list, not a filtered subset.'], ...
        sum(ok)));
    return;
end
genes = genes(ok);
stats = stats(ok);

% ---- 1. Which gene set collection? ----------------------------------
[indx1, species] = gui.i_selgenecollection(parentfig);
if isempty(indx1), return; end
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, species, parentfig);
if isempty(setmatrx) || isempty(setnames) || isempty(setgenes)
    return;
end

% ---- 2. Which test? --------------------------------------------------
% One prompt covers method and direction together; the pairs below are the
% combinations worth offering rather than the full cross product.
choices = ["CAMERA - corrected for gene-gene correlation (recommended)"
    "Wilcoxon rank-sum (limma geneSetTest)"
    "CERNO - sensitive to the very top of the list"
    "PAGE - parametric z-score"
    "GSEA - weighted Kolmogorov-Smirnov (slow)"
    "CAMERA - perturbed in either direction (mixed)"];
methodlist = ["camera", "ranksum", "cerno", "page", "gsea", "camera"];
dirlist = ["both", "both", "both", "both", "both", "mixed"];

if isempty(statlabel)
    prompt = 'Select a competitive gene-set test.';
else
    prompt = sprintf('Select a competitive gene-set test. Ranking by %s.', ...
        statlabel);
end
[indx2, tf2] = gui.myListdlg(parentfig, choices, 'Gene-Set Test', ...
    choices(1), false, true, [430, 260], prompt);
if tf2 ~= 1 || isempty(indx2), return; end
method = methodlist(indx2);
direction = dirlist(indx2);

% ---- 3. Run ----------------------------------------------------------
fw = gui.myWaitbar(parentfig, [], false, ...
    sprintf('Running %s gene-set test...', upper(method)));
try
    T = sc_gsettest(stats, genes, setmatrx, setnames, setgenes, ...
        Method=method, Direction=direction, MaxSize=500);
catch ME
    gui.myWaitbar(parentfig, fw, true);
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    T = table.empty;
    return;
end
gui.myWaitbar(parentfig, fw);

if isempty(T)
    gui.myHelpdlg(parentfig, ...
        'No gene set could be tested against this gene list.', '');
    return;
end

% ---- 4. Show ---------------------------------------------------------
nsig = sum(T.FDR < 0.05);
views(1) = struct('Name', sprintf('All sets (%d)', height(T)), 'Table', T);
views(2) = struct('Name', sprintf('FDR<0.05 (%d)', nsig), ...
    'Table', T(T.FDR < 0.05, :));
gui.TableViewerApp(views, parentfig, defname);
end
