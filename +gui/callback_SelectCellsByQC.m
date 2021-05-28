function callback_SelectCellsByQC(src)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    listitems={'SC_QCFILTER (Basic QC)',...
        'Remove Genes by Expression',...
        'Remove Mt-genes',...
        'Select & remove genes',...
        'Remove Cells by Mt-reads Ratio & Library Size'};
        %    'Remove Cells by Dropout Rate & Expression Mean'};
    [indx,tf] = listdlg('PromptString',{'Select Filter',...
    '',''},'SelectionMode','single',...
    'ListString',listitems,'ListSize',[200 300]);
    if tf~=1, return; end
        switch indx
            case 1
                fw=gui.gui_waitbar;
                sce=sce.qcfilter;
                gui.gui_waitbar(fw);
            case 2
                answer=inputdlg('Expressed in less than % of cells','Remove Genes',[1 40],{'5'});
                if isempty(answer), return; end
                if iscell(answer)
                    a=str2double(answer{1});
                    if a>0 && a<100
                        sce=sce.selectkeepgenes(1,a/100);
                    end
                end
            case 3
                sce=sce.rmmtgenes;
            case 4           % remove selected genes
                gsorted=sort(sce.g);
                [idx]=gui.gui_selmultidlg(gsorted);
                if isempty(idx), return; end
                if isscalar(idx) && idx==0
                    helpdlg('No gene selected.');
                    return;
                else
                    [~,i]=ismember(gsorted(idx),sce.g);
                    answer1 = questdlg('Remove selected genes?');
                    if strcmpi(answer1,'Yes')        
                        sce.g(i)=[];
                        sce.X(i,:)=[];
                    else
                        return;
                    end
                end
            case 5
                i=startsWith(sce.g,'mt-','IgnoreCase',true);
                if ~any(i) 
                    disp('No mt genes found.');
                    return;
                end
                lbsz=sum(sce.X,1);
                lbsz_mt=sum(sce.X(i,:),1);
                cj=lbsz_mt./lbsz;
                if issparse(cj), cj=full(cj); end
                ttxtj="mtDNA%";
                ci=sum(sce.X);
                if issparse(ci), ci=full(ci); end
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
        end
    guidata(FigureHandle,sce);
end
                
%             case 6
%                 i=startsWith(sce.g,'mt-','IgnoreCase',true);
%                 if ~any(i) 
%                     disp('No mt genes found.');
%                     return;
%                 end
%                 rdrop=sum(sce.X==0,2)./size(sce.X,2);
%                 rmean=-log2(mean(sce.X,2)+0.1);
%                 if issparse(rmean), rmean=full(rmean); end
%                 if issparse(rdrop), rdrop=full(rdrop); end
%                 
%                 ttxti="Dropout Rate";
%                 ttxtj="-log2(Expression Mean+0.1)";
%                 a=maxk(rdrop,10);
%                 idx=gui.gui_setranges2(rdrop(:),rmean(:),[0 a(end)],...
%                         [0 max(rmean)],ttxti,ttxtj);
%                 if any(~idx)
%                     answer = questdlg(sprintf('Remove %d cells?',sum(~idx)));
%                     if strcmpi(answer,'Yes')
%                         sce=sce.removecells(~idx);                        
%                     end
%                 end                
