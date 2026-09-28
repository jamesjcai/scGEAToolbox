function callback_GSEAPhenotype(src, ~)
%CALLBACK_GSEAPHENOTYPE  Standard GSEA with sample-label permutation.
%
%   gui.callback_GSEAPhenotype(src)
%
% Runs RUN.ML_GSEA from the GUI: pick the two phenotypes to compare, pick the
% variable that says which biological sample each cell came from, and the
% cells are pooled into one pseudobulk profile per sample before the test.
%
% THE POOLING IS NOT A CONVENIENCE, it is what makes the test valid. This
% is the one enrichment method here whose null comes from relabelling the
% samples, and that null only means anything if the columns are independent
% replicates. Cells are not: thousands of cells from one donor are one
% donor's measurement repeated, and handing them over as columns would
% report nearly every gene set as significant. So the sample variable is
% required, and every sample has to sit entirely within one phenotype -- a
% sample carrying both labels has no label to permute.
%
% GUI.I_RUNGSETTEST is the other door into set testing and takes a ranked
% gene list instead. Use that one when there are no replicates to permute;
% it permutes genes, which needs no replication but assumes genes are
% independent. See RUN.ML_GSEA for what that assumption costs.
%
% See also RUN.ML_GSEA, GUI.I_RUNGSETTEST, SC_GSETTEST.

numPerm = 1000;

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~isempty(FigureHandle) && pkg.i_isvalid(FigureHandle) && FigureHandle.Visible == "on"
    figure(FigureHandle);
    cleanupObj = onCleanup(@() gui.i_raisefig(FigureHandle));
end

% ---- 1. Which two phenotypes? ----------------------------------------
[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end

% ---- 2. Which variable identifies a biological sample? ---------------
[thiss, slabel] = gui.i_select1class(sce, false, ...
    'Which variable identifies a biological sample (donor, animal, well)?', ...
    'Batch ID', FigureHandle);
if isempty(thiss), return; end

sampleid = string(thiss(:));
inuse = i1 | i2;
isPosCell = i2(inuse);                  % group 2 is the positive class
sampleid = sampleid(inuse);

[pools, ~, poolof] = unique(sampleid);
numPool = numel(pools);
poolIsPos = false(numPool, 1);
mixed = strings(0, 1);
for k = 1:numPool
    inpool = poolof == k;
    if all(isPosCell(inpool))
        poolIsPos(k) = true;
    elseif any(isPosCell(inpool))
        mixed(end+1, 1) = pools(k);     %#ok<AGROW> -- a handful at most
    end
end

if ~isempty(mixed)
    gui.myErrordlg(FigureHandle, sprintf( ...
        ['%s of "%s" (%s) holds cells from both "%s" and "%s". A ' ...
        'phenotype permutation test needs one label per sample, so pick ' ...
        'a variable whose levels sit entirely on one side of the ' ...
        'comparison.'], pkg.i_plural(numel(mixed), 'level'), slabel, ...
        strjoin(mixed(1:min(5, numel(mixed))), ', '), string(cL1), ...
        string(cL2)));
    return;
end

n1 = sum(poolIsPos);
n2 = numPool - n1;
if min(n1, n2) < 2
    gui.myErrordlg(FigureHandle, sprintf( ...
        ['"%s" has %d sample(s) and "%s" has %d, by "%s". Label ' ...
        'permutation needs at least two in each group, and realistically ' ...
        'seven. With no replication, use Gene-Set Test on a ranked gene ' ...
        'list instead.'], string(cL2), n1, string(cL1), n2, slabel));
    return;
end

% ---- 3. Which gene set collection? -----------------------------------
[indx1, species] = gui.i_selgenecollection(FigureHandle);
if isempty(indx1), return; end
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(indx1, species, FigureHandle);
if isempty(setmatrx) || isempty(setnames) || isempty(setgenes)
    return;
end

% ---- 4. Pseudobulk, then test ----------------------------------------
fw = gui.myWaitbar(FigureHandle, [], false, 'Pooling cells into samples...');
try
    X = sce.X(:, inuse);
    P = zeros(size(X, 1), numPool);
    for k = 1:numPool
        P(:, k) = full(sum(X(:, poolof == k), 2));
    end
    % Summed counts are on whatever depth each pool happened to have, so
    % they are normalised before being compared across samples. The log is
    % what the ranking metric expects: a difference of means over a spread
    % on the raw scale is dominated by the few most expressed genes.
    P = log1p(sc_norm(P));
    genes = sce.g;
    detected = any(P > 0, 2);
    P = P(detected, :);
    genes = genes(detected);

    gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Permuting %d sample labels...', numPool), 0.2);
    T = run.ml_GSEA(P, genes, poolIsPos, setmatrx, setnames, setgenes, ...
        NumPerm=numPerm, Verbose=false);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(T)
    gui.myHelpdlg(FigureHandle, ...
        'No gene set could be tested against these samples.', '');
    return;
end

% ---- 5. Show ----------------------------------------------------------
% Up and down are split because the sign is the answer here: ES > 0 means
% the set is enriched in the group named in the tab, and a single ranking
% of both mixed together hides which side each hit came from.
nsig = sum(T.FDR < 0.05);
views(1) = struct('Name', sprintf('All sets (%d)', height(T)), 'Table', T);
views(2) = struct('Name', sprintf('FDR<0.05 (%d)', nsig), ...
    'Table', T(T.FDR < 0.05, :));
views(3) = struct('Name', sprintf('Up in %s', string(cL2)), ...
    'Table', T(T.ES > 0, :));
views(4) = struct('Name', sprintf('Up in %s', string(cL1)), ...
    'Table', T(T.ES < 0, :));

defname = sprintf('%s_vs_%s_GSEA', ...
    matlab.lang.makeValidName(string(cL2)), ...
    matlab.lang.makeValidName(string(cL1)));
gui.TableViewerApp(views, FigureHandle, defname);

end
