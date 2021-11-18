function [requirerefresh,highlightindex]=callback_SelectCellsByQC(src)
    requirerefresh=true;
    highlightindex=[];
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    % 'SC_QCFILTER (QC Preserves Lowly-expressed Cells/Genes)',...
    listitems={'SC_QCFILTER (Basic QC for Cells/Genes)',...        
        'Remove Genes by Expression',...
        'Remove Genes by Name',...
        'Remove Mt-genes',...
        '------------------------------------------------',...
        'Library Size vs. Mt-reads Ratio',...
        'Library Size vs. Number of Genes',...
        'Abundant lncRNAs vs. Number of Genes',...
        '------------------------------------------------',...
        'QC Metrics in Violin Plots'};
    
%        '------------- Experimental Options -------------',...
%        'Remove ambient RNA contamination (R required)'};

    [indx,tf] = listdlg('PromptString',{'Select Filter',...
                '',''},'SelectionMode','single',...
                'ListString',listitems,'ListSize',[250 300]);
    if tf~=1, return; end
    switch listitems{indx}
        case 'SC_QCFILTER (Basic QC for Cells/Genes)'   % basic QC
%             answer=questdlg({'Library Size > 1000','mtDNA Ratio < 10%',...
%                                'Gene''s min_cells_nonzero > 5%'});
            
            answer3=questdlg('Relaxed or Strigent?',...
                'Cutoff Settings','Relaxed (keep more cells/genes)',...
                'Strigent (remove more cells/genes)','Relaxed (keep more cells/genes)');
            switch answer3
                case 'Relaxed (keep more cells/genes)'
                    definput = {'500','0.15','0.01'};
                case 'Strigent (remove more cells/genes)'
                    definput = {'1000','0.10','0.05'};
                otherwise
                    requirerefresh=false;
                    return;
            end        
                    
            prompt = {'Library size (e.g., 500 or 1000):','mtDNA ratio (e.g., 15% or 10%):',...
                      'Gene''s min_nonzero_cells (e.g., 0.01 or 0.05, 10 or 50):'};
            dlgtitle = 'QC cutoff';
            dims = [1 55];
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(answer), return; end
            try
                libsize=str2double(answer{1});
                mtratio=str2double(answer{2});
                min_cells_nonzero=str2double(answer{3});
                assert((libsize>100)&&(libsize<100000));
                assert((mtratio>=0) && (mtratio<=1));
                assert((min_cells_nonzero>=0 && min_cells_nonzero<=1)||(min_cells_nonzero>1 && min_cells_nonzero<sce.NumCells));
            catch
                requirerefresh=false;
                errordlg('Invalid input(s).');
                return;
            end            
            [whitelist]=i_selectwhitelist(sce);            
            fw=gui.gui_waitbar;
            sce=sce.qcfilterwhitelist(libsize,mtratio,min_cells_nonzero,whitelist);
            gui.gui_waitbar(fw);
        %   [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
        %   [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);
 
        case 'Remove Genes by Expression'     % remove genes by expression
            answer=inputdlg('Expressed in less than % of cells',...
                'Remove Genes',[1 40],{'5'});
            if isempty(answer), return; end
            if iscell(answer)
                a=str2double(answer{1});
                if a>0 && a<100
                    fw = gui.gui_waitbar; 
                    sce=sce.selectkeepgenes(1,a/100);
                    gui.gui_waitbar(fw);
                end
            end
        case 'Remove Genes by Name'        % remove selected genes
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
            	    fw = gui.gui_waitbar;                    
                    sce.g(i)=[];
                    sce.X(i,:)=[];
                    gui.gui_waitbar(fw);
                else
                    return;
                end
            end
        case 'Remove Mt-genes'        % remove mt-genes
            sce=sce.rmmtgenes;            
        case '------------------------------------------------'
            requirerefresh=false;
            return;
        case 'Library Size vs. Mt-reads Ratio'      % mt-ratio vs. library size
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
        case 'Library Size vs. Number of Genes'
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
        case 'Abundant lncRNAs vs. Number of Genes'    % 'Abundant lncRNAs vs. Number of Genes'
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
        case 'QC Metrics in Violin Plots'   % view QC metrics violin
            gui.sc_qcviolin(sce.X,sce.g);
            requirerefresh=false;
            return;         
%         case 11
%             fw = gui.gui_waitbar;
%             [Xdecon,contamination]=run.decontX(sce);
%             sce.X=Xdecon;
%             guidata(FigureHandle,sce);
%             gui.gui_waitbar(fw);
%             figure;
%             gui.i_stemscatter(sce.s,contamination);
%             zlabel('Contamination rate')
%             title('Ambient RNA contamination')
%             pause(1)
%             helpdlg('Contamination removed.')
%             requirerefresh=false;
%             return;
        otherwise
            requirerefresh=false;
            return;
    end
    
    if ismember(indx,[6 7 8])
        if any(~idx)
            answer = questdlg(sprintf('Remove or highlight %d cells?',sum(~idx)),...
                '','Remove','Highlight','Cancel','Remove');
            switch answer
                case 'Remove'
                    sce=sce.removecells(~idx);
                case 'Highlight'              
                    highlightindex=zeros(1,length(idx));
                    highlightindex(~idx)=1;
                    requirerefresh=false;
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
        end
    end
    guidata(FigureHandle,sce);
end
                


function [whitelist]=i_selectwhitelist(sce)
    whitelist=[];        
    answer = questdlg('Genes in whitelist will not be removed. Select whitelist genes?',...
                'Whitelist Genes','Yes','No','Cancel','No');
    switch answer
        case 'Yes'
            [gsorted]=gui.i_sortgenenames(sce);
            if isempty(gsorted), return; end
            [idx]=gui.i_selmultidlg(gsorted);
            if isempty(idx), return; end
            if isscalar(idx) && idx==0, return; end
            whitelist=gsorted(idx);
        case 'No'
            return;
        case 'Cancel'
            return;
        otherwise
            return;
    end
end