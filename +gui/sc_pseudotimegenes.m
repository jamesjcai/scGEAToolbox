function sc_pseudotimegenes(sce, t, parentfig)

if nargin<3, parentfig=[]; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
t = t(:);

[K, genemode] = gui.i_gethvgnum(sce, parentfig);
% Dismissing the dialog used to leave K empty and USEHVGS false, which fell
% through to the all-genes branch and ran the analysis anyway.
if isempty(K), return; end

if genemode == "all"
    sceX = sce.X;
    sceg = sce.g;
else
    T = sc_analyticfit(sce.X, sce.g);
    % HEIGHT(T) rather than SCE.NUMGENES: the ranker drops genes that
    % are zero in every cell, so its table can be shorter than the SCE
    % and indexing out to NUMGENES throws.
    glist = T.genes(1:min([K, height(T)]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y), error('Runtime error.'); end
    sceX = sce.X(idx, :);
    sceg = sce.g(idx);
    if genemode == "hvg+markers"
        [sceX, sceg] = pkg.i_appendgenes(sceX, sceg, ...
            pkg.i_getmarkerwhitelist(sce.g, sce.X), sce.X, sce.g);
    end
end

X = sceX;
try
    if issparse(X), X = full(X); end
catch ME
    disp(ME.message);
    disp('Keep using sparse X.');
end


% Four choices: gui.myQuestdlg cannot carry more than two plus Cancel.
methods = {'Spearman Correlation', 'Distance Correlation', ...
    'MIC (MINE)', 'MAS (MINE)'};
prompt = ['Rank genes against pseudotime. Spearman finds monotone ' ...
    'trends; distance correlation and MIC find any dependence; MAS ' ...
    'finds non-monotonic ones, such as genes that rise and then fall.'];
[sel, ok] = gui.myListdlg(parentfig, methods, 'Ranking Method', ...
    methods{1}, false, true, [420, 170], prompt);
if ~ok || isempty(sel), return; end
answer = methods{sel};

switch answer
    case 'Spearman Correlation'
        fw = gui.myWaitbar(parentfig);
        r = corr(t, X.', 'type', 'spearman');
        gui.myWaitbar(parentfig, fw);
    case 'Distance Correlation'
        Xn = log1p(sc_norm(X));
        r = i_scoregenes(parentfig, Xn, t, @pkg.e_distcorr);
    case {'MIC (MINE)', 'MAS (MINE)'}
        Xn = log1p(sc_norm(X));
        % The grid search grows faster than linearly in the number of
        % cells: on normalised data about 0.04 s a gene at 1000 cells and
        % 0.13 s at 2000. Fewer cells cost power, though -- of 20 weak
        % transient genes among 200 noise genes, MAS put 15 in its top 20
        % at 1000 cells and 19 at 2000. Cells evenly spaced along
        % pseudotime keep its whole range.
        maxMineCells = 2000;
        [~, order] = sort(t);
        keep = order(unique(round(linspace(1, numel(t), ...
            min(numel(t), maxMineCells)))));
        if answer == "MIC (MINE)"
            field = "mic";
        else
            field = "mas";
        end
        r = i_scoregenes(parentfig, Xn(:, keep), t(keep), ...
            @(x, y) i_minestat(x, y, field));
    otherwise
        % myListdlg returns only items from METHODS.
        return;
end

% By magnitude: a gene falling along pseudotime (Spearman r near -1) is as
% much a trend gene as a rising one. The other scores are non-negative.
[~, idxp] = maxk(abs(r), 10);
selectedg = sceg(idxp);
try
        hx=gui.myFigure(parentfig);
        hFig=hx.FigHandle;
        hFig.Position(3) = hFig.Position(3) * 1.8;
        hx.show(parentfig);
        figure(hFig);
        gui.i_plot_pseudotimeseries(log1p(X), ...
            sceg, t, selectedg);
        title(sprintf('Top %d genes by %s', numel(selectedg), answer));
    catch ME
        if exist('psf1', 'var') && ishandle(hFig)
            close(hFig);
        end
        gui.myErrordlg(parentfig, ME.message,'','modal');
end

end

function r = i_scoregenes(parentfig, Xn, t, scorefun)
% Score each gene (row of XN) against pseudotime T with SCOREFUN(x, t).
% Genes are split across an already-open parallel pool, if there is one;
% none is started here, since that can take longer than the scoring. The
% chunks exist so the waitbar can move between parfor loops.
numGenes = size(Xn, 1);
t = t(:);
numWorkers = i_poolsize();
chunkSize = max(100, 25*numWorkers);
r = zeros(numGenes, 1);
fw = gui.myWaitbar(parentfig);
for first = 1:chunkSize:numGenes
    idx = first:min(first + chunkSize - 1, numGenes);
    Xc = Xn(idx, :);
    rc = zeros(numel(idx), 1);
    parfor (k = 1:numel(idx), numWorkers)
        rc(k) = scorefun(Xc(k, :).', t);
    end
    r(idx) = rc;
    gui.myWaitbar(parentfig, fw, false, '', '', idx(end)/numGenes);
end
gui.myWaitbar(parentfig, fw);
end

function n = i_poolsize()
% Workers in the current pool, or 0 -- which makes parfor run serially and
% never auto-create a pool. GCP is absent without Parallel Computing
% Toolbox, so its failure means no pool.
n = 0;
try
    p = gcp("nocreate");
    if ~isempty(p)
        n = p.NumWorkers;
    end
catch
    % No Parallel Computing Toolbox: run serially.
end
end

function v = i_minestat(x, t, field)
% One MINE statistic. Noise genes score MIC and MAS of about 0.1, not 0:
% the grid search always finds some structure, so read ranks, not values.
s = run.ml_mine(x, t);
v = s.(field);
end
