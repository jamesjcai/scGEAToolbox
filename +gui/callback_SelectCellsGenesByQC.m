function callback_SelectCellsGenesByQC(src)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    listitems={'SC_QCFILTER (General QC)','Remove Genes by Expression',...
        'Remove Mt-genes',...
        'Remove Cells by Mt-reads Ratio & Library Size'};
    [indx,tf] = listdlg('PromptString',{'Select Filter',...
    '',''},'SelectionMode','single',...
    'ListString',listitems,'ListSize',[200 300]);
    if tf~=1, return; end
        switch indx
            case 1
                sce=sce.qcfilter;
            case 2
                answer=inputdlg('Expressed in less than % of cells','Remove Genes',[1 40],{'5'});
                if isempty(answer), return; end
                if iscell(answer)
                    a=str2double(answer{1});
                    if a>0 && a<100
                        sce=sce.selectgenes(1,a/100);
                    end
                end
            case 3
                sce=sce.rmmtgenes;
            case 4
                i=startsWith(sce.g,'mt-','IgnoreCase',true);
                if ~any(i) 
                    disp('No mt genes found.');
                    return;
                end
                lbsz=sum(sce.X,1);
                lbsz_mt=sum(sce.X(i,:),1);
                cj=lbsz_mt./lbsz;
                ttxtj="mtDNA%";
                ci=sum(sce.X);
                ttxti="Library Size";
                a=maxk(ci,10);
                idx=gui.gui_setranges2(ci',cj',[0 a(end)],...
                        [0 0.1],ttxti,ttxtj);
                if any(~idx)
                    answer = questdlg(sprintf('Remove %d cells?',sum(~idx)));
                    if strcmpi(answer,'Yes')
                        sce=sce.removecells(~idx);                        
                    end
                end
            otherwise
                
        end
    guidata(FigureHandle,sce);
end
