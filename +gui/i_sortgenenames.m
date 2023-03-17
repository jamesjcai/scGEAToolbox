function [gsorted]=i_sortgenenames(sce)
        gsorted=[];
        answer2 = questdlg('How to sort genes?','Sort Genes',...
            'Alphabetic','Average Expression','% of Nonzero Cells','Alphabetic');
        switch answer2
            case 'Alphabetic'
                gsorted=natsort(sce.g);
            case 'Average Expression'
                [T]=sc_genestats(sce.X,sce.g);
                [~,idx]=sort(T.Mean,'descend');
                gsorted=sce.g(idx);                
            case '% of Nonzero Cells'
                [T]=sc_genestats(sce.X,sce.g);
                [~,idx]=sort(T.Dropout_rate,'ascend');
                gsorted=sce.g(idx);
            otherwise
                return;
        end
end