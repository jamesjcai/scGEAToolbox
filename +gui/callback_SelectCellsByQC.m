function [requirerefresh,highlightindex]=callback_SelectCellsByQC(src)
    requirerefresh=true;
    highlightindex=[];
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    listitems={'SC_QCFILTER (Basic QC)',...
        'QC metrics in violin plot',...
        'Remove Genes by Expression',...
        'Remove Mt-genes',...
        'Select & remove genes',... 
        '------------------------------------------------',...
        'Library Size vs. Mt-reads Ratio',...
        'Library Size vs. Number of Genes',...
        'Reads in Abundant lncRNAs vs. Number of Genes'};
%        '------------- Experimental Options -------------',...
%        'Remove ambient RNA contamination (R required)'};

    [indx,tf] = listdlg('PromptString',{'Select Filter',...
                '',''},'SelectionMode','single',...
                'ListString',listitems,'ListSize',[250 300]);
    if tf~=1, return; end
    switch indx
        case 1   % basic QC
            fw=gui.gui_waitbar;            
            sce=sce.qcfilter;
            gui.gui_waitbar(fw);
        case 2   % view QC metrics violin
            gui.sc_qcviolin(sce.X,sce.g);
            requirerefresh=false;
            return;            
        case 3     % remove genes
            answer=inputdlg('Expressed in less than % of cells',...
                'Remove Genes',[1 40],{'5'});
            if isempty(answer), return; end
            if iscell(answer)
                a=str2double(answer{1});
                if a>0 && a<100
                    sce=sce.selectkeepgenes(1,a/100);
                end
            end
        case 4        % remove mt-genes
            sce=sce.rmmtgenes;
        case 5          % remove selected genes
            gsorted=sort(sce.g);
            [idx]=gui.i_selmultidlg(gsorted);
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
        case 6
            % -----------
            requirerefresh=false;
            return;
        case 7      % mt-ratio vs. library size
            i=startsWith(sce.g,'mt-','IgnoreCase',true);
            if ~any(i) 
                disp('No mt genes found.');
                return;
            end
            lbsz=sum(sce.X,1);
            lbsz_mt=sum(sce.X(i,:),1);
            cj=100*(lbsz_mt./lbsz);
            if issparse(cj), cj=full(cj); end
            ttxtj="mtDNA%";

            ci=sum(sce.X);
            if issparse(ci), ci=full(ci); end
            ttxti="Library Size";
            a=maxk(ci,10);                
            idx=gui.i_setranges2(ci',cj',[0 a(end)],...
                    [0 15],ttxti,ttxtj);
        case 8
            cj=sum(sce.X>0,1);
            if issparse(cj), cj=full(cj); end
            ttxtj="Number of Detected Genes";
            ci=sum(sce.X,1);
            if issparse(ci), ci=full(ci); end
            ttxti="Library Size";
            a=maxk(ci,10);
            b=maxk(cj,10);
            idx=gui.i_setranges2(ci',cj',[0 a(end)],...
                    [0 b(end)],ttxti,ttxtj);
        case 9
            % remove cells with a high fraction of nuclear lncRNA transcripts 
            % (Malat1, Meg3 and Kcnq10t1)
            % https://www.frontiersin.org/articles/10.3389/fncel.2020.00065/full#h3
            cj=sum(sce.X>0,1);
            if issparse(cj), cj=full(cj); end
            ttxtj="Number of Detected Genes";
                       
            idx=matches(upper(sce.g),upper({'Malat1', 'Meg3', 'Kcnq1ot1'}));
            if ~any(idx) 
                disp('{Malat1,Meg3,Kcnq1ot1} not found.');
                return;
            end
            ci=sum(sce.X(idx,:),1);            
            if issparse(ci), ci=full(ci); end
            ttxti="Reads in lncRNAs (Malat1,Meg3,Kcnq1ot1)";
            
            a=maxk(ci,10);
            b=maxk(cj,10);
            idx=gui.i_setranges2(ci',cj',[0 a(end)],...
                    [0 b(end)],ttxti,ttxtj);
        case 10
            % ----------
            requirerefresh=false;
            return;
        case 11
            fw = gui.gui_waitbar;
            [Xdecon,contamination]=run.decontX(sce);
            sce.X=Xdecon;
            guidata(FigureHandle,sce);
            gui.gui_waitbar(fw);
            figure;
            gui.i_stemscatter(sce.s,contamination);
            zlabel('Contamination rate')
            title('Ambient RNA contamination')
            pause(1)
            helpdlg('Contamination removed.')
            requirerefresh=false;
            return;
    end
    if ismember(indx,[7 8])
        if any(~idx)
            answer = questdlg(sprintf('Remove or highlight %d cells?',sum(~idx)),...
                '','Remove','Highlight','Cancel','Remove');
            switch answer
                case 'Remove'
                    sce=sce.removecells(~idx);
                case 'Highlight'              
                    highlightindex=zeros(1,length(idx));
                    highlightindex(~idx)=1;
                case 'Cancel'
                    return;
                otherwise
                    return;
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
%                 idx=gui.i_setranges2(rdrop(:),rmean(:),[0 a(end)],...
%                         [0 max(rmean)],ttxti,ttxtj);
%                 if any(~idx)
%                     answer = questdlg(sprintf('Remove %d cells?',sum(~idx)));
%                     if strcmpi(answer,'Yes')
%                         sce=sce.removecells(~idx);                        
%                     end
%                 end                
