function callback_CellHeatMap(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);


% species=questdlg('Which species?','Select Species','Mouse','Human','Mouse');
% switch lower(species)
%     case 'human'
%         stag='hs';
%     case 'mouse'
%         stag='mm';
%     otherwise
%         return;
% end



hFigure=figure;
        UitoolbarHandle = uitoolbar('Parent', hFigure);

pkg.i_addbutton2fig(UitoolbarHandle)
[~,idx]=sort(sce.c);
X=sce.X(:,idx);

imagesc(log10(1+log10(1+X(1:1000,:))));
xlabel("Genes")
ylabel("Cells")


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
        [c,cL]=grp2idx(thisc);        
        sce.c=c;
    end

end
