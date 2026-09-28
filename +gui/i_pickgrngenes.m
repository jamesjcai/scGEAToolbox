function [glist, useAll] = i_pickgrngenes(genes, parentfig)
%I_PICKGRNGENES Choose the genes a gene network covers.
%   [GLIST, USEALL] = gui.i_pickgrngenes(GENES, parentfig) asks whether to
%   select genes from a list, paste gene names, or use all genes. GLIST is
%   the chosen genes (at least two) as spelled in GENES; USEALL is true
%   when the user chose All Genes, and then GLIST is empty. Both are empty
%   / false if the user cancels.
%
%   See also gui.callback_BuildGeneNetwork, gui.callback_CompareGeneNetwork.

glist = [];
useAll = false;
answer = gui.myQuestdlg(parentfig, sprintf(['Which genes should be ', ...
    'included?\n\nSelect or paste a small set of genes, or use all %d ', ...
    'genes.'], numel(genes)), 'Genes', ...
    {'Select Genes...', 'Paste Genes...', 'All Genes'}, 'Select Genes...');
switch answer
    case 'Select Genes...'
        gsorted = natsort(genes);
        idx = gui.i_selmultidialog(gsorted, [], parentfig);
        if isempty(idx) || isequal(idx, 0), return; end
        glist = gsorted(idx);
    case 'Paste Genes...'
        glist = i_pastegenes(genes, parentfig);
    case 'All Genes'
        useAll = true;
        return;
    otherwise
        return;
end
if isempty(glist), return; end
if numel(glist) < 2
    gui.myWarndlg(parentfig, ['A network needs at least two genes. ', ...
        'Select more genes and try again.']);
    glist = [];
end
end


function glist = i_pastegenes(genes, parentfig)
% Genes pasted as names, matched case-insensitively to GENES, or [] if
% cancelled. Names that match no gene are reported and dropped.
glist = [];
tg = gui.i_inputgenelist(strings(0, 1), false, parentfig);
if isempty(tg), return; end
[y, ix] = ismember(upper(tg), upper(genes));
glist = genes(ix(y));
nMissing = nnz(~y);
if nMissing > 0
    gui.myWarndlg(parentfig, sprintf('%s not found: %s', ...
        pkg.i_plural(nMissing, 'gene'), strjoin(tg(~y), ', ')));
end
end
