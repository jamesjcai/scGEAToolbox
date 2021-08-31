function [gsorted]=i_sortgenenames(sce)
        gsorted=[];
        answer2 = questdlg('How to sort gene names?','Sort by',...
            'Alphabetic','Expression Mean','Dropoff Rate','Alphabetic');
        switch answer2
            case 'Alphabetic'
                gsorted=sort(sce.g);
            case 'Expression Mean'
                [T]=sc_genestats(sce.X,sce.g);
                [~,idx]=sort(T.Dropout_rate);
                gsorted=sce.g(idx);                
            case 'Dropoff Rate'
                [T]=sc_genestats(sce.X,sce.g);
                [~,idx]=sort(T.Dropout_rate);
                gsorted=sce.g(idx);
            otherwise
                return;
        end
end