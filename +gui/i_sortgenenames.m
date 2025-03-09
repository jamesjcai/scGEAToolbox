function [gsorted] = i_sortgenenames(sce, parentfig)
if nargin < 2, parentfig = []; end
gsorted = [];
answer2 = gui.myQuestdlg(parentfig, 'How to sort genes?', 'Sort Genes', ...
    {'Alphabetic', 'Average Expression'}, 'Alphabetic');

switch answer2
    case 'Alphabetic'
        gsorted = natsort(sce.g);
    case 'Average Expression'
        %X=sce.X;
        %if issparse(X)
        %    try                
        %        X=full(X);
        %    catch
        %    end
        %end
        % [T] = sc_genestats(sce.X, sce.g);
        [~, idx] = sort(mean(sce.X,2), 'descend');
        gsorted = sce.g(idx);
    case '% of Nonzero Cells'
        %[T] = sc_genestats(sce.X, sce.g);
        %[~, idx] = sort(T.Dropout_rate, 'ascend');
        %gsorted = sce.g(idx);
        tic;
        X=sce.X;
        if issparse(X)
           try                
               X=full(X);
           catch
           end
        end
        %tic;        
            [~, idx] = sort(mean(X,2), 'descend');
            gsorted = sce.g(idx);
            toc;
            X=X(idx,:);
            [~, idx] = sort(sum(X>0,2), 'descend');
            gsorted = gsorted(idx);
        %toc;
    otherwise
        return;
end
end