function [glist] = i_selectngenes(sce, predefinedlist, parentfig)



if nargin < 2, predefinedlist = []; end
if nargin < 3, parentfig = []; end

if ~isempty(parentfig)
    figure(parentfig);
    cleanupObj = onCleanup(@() figure(parentfig));
end

% internal function used by callback_BuildGeneNetwork
glist = [];
if isa(sce, 'SingleCellExperiment')
    gsorted = natsort(sce.g);
elseif isa(sce, 'SingleCellExperiment2')
    gsorted = natsort(sce.geneAnn.names);
elseif isstring(sce)
    genelist = sce;
    gsorted = natsort(genelist);
end

if ~isempty(predefinedlist)
    predefinedlist = gsorted(matches(gsorted, predefinedlist, ...
        'IgnoreCase', true));
end

answer = gui.myQuestdlg(parentfig, 'Select genes from list or paste gene names?', ...
    'Select/Paste Genes', {'Select', 'Paste', 'Cancel'}, 'Select');
switch answer
    case 'Select'
        if isa(sce, 'SingleCellExperiment')
            [gsorted] = gui.i_sortgenenames(sce, parentfig);
            if isempty(gsorted), return; end
        elseif isa(sce, 'SingleCellExperiment2')
            [gsorted] = gui.i_sortgenenames(sce, parentfig);
            if isempty(gsorted), return; end
        end

        [idx] = gui.i_selmultidialog(gsorted, predefinedlist, parentfig);
        if isempty(idx), return; end
        if idx == 0, return; end
        glist = gsorted(idx);
    case 'Paste'
        rng("shuffle");
        n = length(gsorted);
        if isempty(predefinedlist)
            predefinedlist = gsorted(randperm(n, ...
                min([15, n])));
        end
        tg = gui.i_inputgenelist(predefinedlist, false, parentfig);
        if length(tg) >= 1
            [y, ix] = ismember(upper(tg), upper(gsorted));
            glist = gsorted(ix(y));
            a = length(glist) - length(tg);
            if a ~= 0
                if a == 1
                    gui.myWarndlg(parentfig, sprintf('%d gene is not found.', a));
                elseif a > 1
                    gui.myWarndlg(parentfig, sprintf('%d genes are not found.', a));
                end
            end
        end
    case 'Cancel'
        return;
    otherwise
        return;
end

end
