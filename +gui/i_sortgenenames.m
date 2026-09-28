function [gsorted] = i_sortgenenames(sce, parentfig)
if nargin < 2, parentfig = []; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
gsorted = [];
answer2 = gui.myQuestdlg(parentfig, 'How to sort genes?', 'Sort Genes', ...
{'Alphabetic', 'Average Expression', 'Unsorted'}, 'Alphabetic');
if isempty(answer2), return; end

if isa(sce, 'SingleCellExperiment')

    switch answer2
        case 'Alphabetic'
            gsorted = natsort(sce.g);
        case 'Average Expression'
            [~, idx] = sort(mean(sce.X,2), 'descend');
            gsorted = sce.g(idx);
        case 'Unsorted'
            gsorted = sce.g;
        case '% of Nonzero Cells'
            tic;
            X=sce.X;
            if issparse(X)
               try
                   X=full(X);
               catch
                   % keep X sparse if it's too large to densify
               end
            end
                [~, idx] = sort(mean(X,2), 'descend');
                gsorted = sce.g(idx);
                toc;
                X=X(idx,:);
                [~, idx] = sort(sum(X>0,2), 'descend');
                gsorted = gsorted(idx);
        otherwise
            return;
    end

end

end
