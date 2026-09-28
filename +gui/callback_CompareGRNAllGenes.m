function callback_CompareGRNAllGenes(src, i1, i2, name1, name2)
% CALLBACK_COMPAREGRNALLGENES  Compare genome-wide GRNs of two cell groups.
%
% The all-genes branch of Network > Build && Compare Two GRNs...:
% gui.callback_CompareGeneNetwork has already picked the two groups
% (logical I1 and I2, one element per cell) and their names. Here: pick the
% method, build one network per group, align the two by manifold
% alignment and rank the genes whose regulation differs (scTenifoldNet).
% Afterwards: a Q-Q plot of the ranking, then save, optional GSEA, and
% export of the table.
%
% Methods
%   PCR           net.pcrnet on each group, as a single network.
%   Denoised PCR  ten.sctenifoldnet - PCR networks on 10 subsamples of 500
%                 cells per group, merged by tensor decomposition (the full
%                 scTenifoldNet construction; needs Tensor Toolbox). It drops
%                 ribosomal genes, so its networks and table are over fewer
%                 genes than sce.g.
%
% See also gui.callback_CompareGeneNetwork, ten.sctenifoldnet, ten.i_ma,
%   ten.i_dr, gui.callback_scTenifoldNetView.

[FigureHandle, sce] = gui.gui_getfigsce(src);
name1 = string(name1);
name2 = string(name2);

% --- Method
numSubsamples = 10;
subsampleSize = 500;
methods = {'PCR (fast)', ...
    'Denoised PCR (10 subsamples + tensor decomposition; slow)'};
prompt = sprintf(['How should the two networks over all %d genes be ', ...
    'built?\n\n%s: %d cells.  %s: %d cells.\n\nPCR builds one network per ', ...
    'group. Denoised PCR builds PCR networks on %d subsamples of %d ', ...
    'cells per group and merges them; it can take hours and needs Tensor ', ...
    'Toolbox. Either way, the networks are aligned and genes ranked by ', ...
    'how much their regulation differs. [scTenifoldNet, PMID:33336197]'], ...
    numel(sce.g), name1, nnz(i1), name2, nnz(i2), numSubsamples, subsampleSize);
[indx, tf] = gui.myListdlg(FigureHandle, methods, 'Network Method', ...
    1, false, true, [440, 130], prompt);
if tf ~= 1 || isempty(indx), return; end

if indx == 2
    if min(nnz(i1), nnz(i2)) <= subsampleSize
        gui.myWarndlg(FigureHandle, sprintf(['Denoised PCR subsamples ', ...
            '%d cells from each group, but the smaller group has only %d. ', ...
            'Use PCR, or choose larger groups.'], subsampleSize, ...
            min(nnz(i1), nnz(i2))));
        return;
    end
    if ~i_checktensortoolbox(FigureHandle), return; end
end

% --- Build, align, rank
g = sce.g;
fw = gui.myWaitbar(FigureHandle);
try
    switch indx
        case 1
            X = log1p(sc_norm(sce.X));
            X0 = X(:, i1);
            X1 = X(:, i2);
            disp('Constructing networks (1/2) ...')
            A0 = net.pcrnet(X0, 3, false, true, false, false, pkg.i_usegpu(X0));
            disp('Constructing networks (2/2) ...')
            A1 = net.pcrnet(X1, 3, false, true, false, false, pkg.i_usegpu(X1));
            disp('Manifold alignment...')
            [aln0, aln1] = ten.i_ma(0.5*(A0 + A0'), 0.5*(A1 + A1'));
            disp('Differential regulation (DR) detection...')
            T = ten.i_dr(aln0, aln1, g);
            tag = 'pcr';
        case 2
            disp('[T, A0, A1, g] = ten.sctenifoldnet(sce.X(:,i1), sce.X(:,i2), sce.g, ''savegrn'', false);')
            % G is replaced by the genes A0 and A1 are over (no ribosomal genes)
            [T, A0, A1, g] = ten.sctenifoldnet(sce.X(:, i1), sce.X(:, i2), g, ...
                'nsubsmpl', numSubsamples, 'csubsmpl', subsampleSize, ...
                'savegrn', false);
            tag = 'denoisedpcr';
        otherwise
            % myListdlg returns an index into METHODS; nothing else can occur
            error('gui:callback_CompareGRNAllGenes:BadMethod', 'Unknown method.');
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

qqFig = ten.e_mkqqplot(T);
i_saveresult(T, A0, A1, g, name1, name2, tag, FigureHandle, qqFig);
end


function ok = i_checktensortoolbox(FigureHandle)
ok = true;
try
    ten.check_tensor_toolbox;
catch
    gui.i_installtensortoolbox(FigureHandle);
    try
        ten.check_tensor_toolbox;
    catch ME
        gui.myErrordlg(FigureHandle, ME.message);
        ok = false;
    end
end
end


function i_saveresult(T, A0, A1, glist, name1, name2, tag, FigureHandle, qqFig)
% Offer to save the result, run GSEA on the ranking, and export the
% tables. The .mat file holds the DR table T, both networks A0 and A1,
% their gene list glist, and the two group names.
groups = [name1; name2];
defaultFolder = getpref('scgeatoolbox', 'externalwrkpath', '');
if ~isfolder(defaultFolder), defaultFolder = pwd; end
defaultName = matlab.lang.makeValidName(sprintf('tenifoldnet_%s_%s_vs_%s', ...
    tag, name1, name2)) + "_" + string(datetime("now", Format="yyyyMMdd_HHmmss")) + ".mat";
[file, folder] = uiputfile({'*.mat', 'MAT-files (*.mat)'}, ...
    'Save Comparison As', fullfile(defaultFolder, defaultName));
if pkg.i_isvalid(FigureHandle), figure(FigureHandle); end

gseaPrompt = sprintf(['Genes are ranked by how much their regulation ', ...
    'differs between %s and %s.\n\nRun gene set enrichment analysis ', ...
    '(GSEA) on the ranking?'], name1, name2);
if ~isequal(file, 0)
    f1 = fullfile(folder, file);
    save(f1, 'T', 'A0', 'A1', 'glist', 'groups', '-v7.3');
    fprintf('The result has been saved in %s\n', f1);
    gseaPrompt = sprintf('The result has been saved in %s.\n\n%s', f1, gseaPrompt);
end

Tr = [];
if strcmp('Yes', gui.myQuestdlg(FigureHandle, gseaPrompt, 'GSEA', ...
        {'Yes', 'No'}, 'No'))
    % GSEA runs permutations and can take minutes; without a bar the app
    % looked frozen. Closed before any error dialog.
    fw = gui.myWaitbar(FigureHandle, [], [], 'Running GSEA...');
    try
        Tr = ten.e_fgsearun(T);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        Tr = [];
    end
    gui.myWaitbar(FigureHandle, fw);
    if ~isempty(Tr) && strcmp('Yes', gui.myQuestdlg(FigureHandle, ...
            'Group the enriched gene sets into a network of related terms?', ...
            'GSEA', {'Yes', 'No'}, 'Yes'))
        % Into the Q-Q plot's window, with a Back button to it, rather
        % than a third window over the app and the plot.
        gui.myFigure.drawInto(qqFig, @() ten.e_fgseanet(Tr));
    end
end

gui.i_exporttable(T, true, 'Ttenifldnet', 'TenifldNetTable', [], [], FigureHandle);
if ~isempty(Tr)
    gui.i_exporttable(Tr, true, 'Tgseaoutput', 'GSEAResultTable', [], [], FigureHandle);
end
end
