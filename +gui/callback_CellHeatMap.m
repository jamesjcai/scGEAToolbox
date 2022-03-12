function callback_CellHeatMap(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    hFigure=figure;
    UitoolbarHandle = uitoolbar('Parent', hFigure);
    pkg.i_addbutton2fig(UitoolbarHandle,'off',@i_changec,'HDF_object01.gif','Sort cells...');
    pkg.i_addbutton2fig(UitoolbarHandle,'off',@i_changeg,'HDF_object02.gif','Sort genes...');
    
    X=sce.X;
    g=sce.g;
    h=[];
    i_redrawh;
    
    mfolder = fileparts(mfilename('fullpath'));

    function i_changeg(~,~)
        answer=questdlg('Sort genes by?','','Chromosomal Position','Others','Chromosomal Position');
        switch answer
            case 'Chromosomal Position'

                species=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
                switch lower(species)
                    case 'human'
                        stag='genelist_human.txt';
                    case 'mouse'
                        stag='genelist_mouse.txt';
                    otherwise
                        return;
                end

                warning off
                T=readtable(fullfile(mfolder,'..', 'resources', stag));
                warning on
                %c=T.Chromosome_scaffoldName;
                [y,idx]=ismember(upper(g),upper(string(T.GeneName)));
                [~,idx]=sort(idx(y));
            case 'Others'
                [gsorted]=gui.i_sortgenenames(sce);
                [~,idx]=ismember(gsorted,g);
        end
        g=g(idx);
        X=X(idx,:);
        % size(X)
        sce.X=X;
        sce.g=g;
        i_redrawh;
    end

    function i_changec(~,~)
        [thisc,clable,~,newpickclable]=gui.i_select1state(sce);
        if strcmp(clable,'Cell Cycle Phase')
            if length(unique(thisc))>1
                sce.c_cell_cycle_tx=thisc;
            end
        end
        if isempty(thisc), return; end
            if strcmp(clable,'Customized C...')
                clable=gui.i_renamec(clable,sce,newpickclable);
                sce.list_cell_attributes=[sce.list_cell_attributes,clable];
                sce.list_cell_attributes=[sce.list_cell_attributes,thisc];
            end
            [c,~]=grp2idx(thisc);            
            [~,idx]=sort(c);
            X=sce.X;
            X=X(:,idx);
            sce=sce.sortcells(idx);
            i_redrawh;
    end

    function i_redrawh
        if exist('h','var') && ~isempty(h)
            delete(h);
        end
        h=imagesc(10*(log10(1+10*log10(1+X))));
        ylabel("Genes");
        xlabel("Cells");
        title(sce.title)
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcnx};

    end

    function [txt] = i_myupdatefcnx(~, event_obj)
        % pos = event_obj.Position;
        idx = event_obj.DataIndex;
        [row,col] = ind2sub([sce.NumGenes sce.NumCells],idx);
        txt = g(row);      
    end
end  