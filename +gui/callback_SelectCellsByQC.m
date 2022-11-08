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
        'Remove Ribosomal Genes',...
        '------------------------------------------------',...
        'Library Size vs. Mt-reads Ratio',...
        'Library Size vs. Number of Genes',...
        'Abundant lncRNAs vs. Number of Genes',...
        '------------------------------------------------',...
        'QC Metrics in Violin Plots'};
    
%        '------------- Experimental Options -------------',...
%        'Remove ambient RNA contamination (R required)'};

    [indx,tf] = listdlg('PromptString',{'Select Filter'},'SelectionMode','single',...
                'ListString',listitems,'ListSize',[250 300]);
    if tf~=1
        requirerefresh=false;
        return;
    end
    switch listitems{indx}
        case 'SC_QCFILTER (Basic QC for Cells/Genes)'   % basic QC
%             answer=questdlg({'Library Size > 1000','mtDNA Ratio < 10%',...
%                                'Gene''s min_cells_nonzero > 5%'});
            
            answer3=questdlg('Relaxed or Strigent?',...
                'Cutoff Settings','Relaxed (keep more cells/genes)',...
                'Strigent (remove more cells/genes)','Relaxed (keep more cells/genes)');
            switch answer3
                case 'Relaxed (keep more cells/genes)'
                    definput = {'500','0.20','10','200'};
                case 'Strigent (remove more cells/genes)'
                    definput = {'1000','0.15','15','500'};
                otherwise
                    requirerefresh=false;
                    return;
            end        
                    
            prompt = {'Library size (e.g., 500 or 1000):',...
                      'mtDNA ratio (e.g., 0.15=15% or 0.10=10%):',...
                      'Gene''s min_nonzero_cells (e.g., 0.01 or 0.05, 10 or 50):',...
                      'Number of genes (e.g., 200 or 500):'};
            dlgtitle = 'QC Cutoffs';
            dims = [1 65];
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(answer)
                requirerefresh=false;
                return; 
            end
            try
                libsize=str2double(answer{1});
                mtratio=str2double(answer{2});
                min_cells_nonzero=str2double(answer{3});
                numgenes=str2double(answer{4});
                assert((libsize>0)&&(libsize<intmax));                
                assert((mtratio>=0.0) && (mtratio<=1.0));
                assert((min_cells_nonzero>=0 && min_cells_nonzero<=1)||(min_cells_nonzero>1 && min_cells_nonzero<sce.NumCells));
                assert((numgenes>0)&&(numgenes<intmax));
            catch
                requirerefresh=false;
                errordlg('Invalid input(s).');
                return;
            end

            [whitelist]=gui.i_selectwhitelist(sce);
            %[whitelist]=gui.i_selectngenes(sce);

            if isnumeric(whitelist) 
                if whitelist==0
                    requirerefresh=false;
                    return;
                end
            end
            % when isempty(whitelist), continue...

            fw=gui.gui_waitbar;
            sce=sce.qcfilterwhitelist(libsize,mtratio,...
                min_cells_nonzero,numgenes,whitelist);
            gui.gui_waitbar(fw);
        %   [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
        %   [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);
 
        case 'Remove Genes by Expression'     % remove genes by expression
            answer=inputdlg('Expressed in less than % of cells (e.g., 0.05=5%) or # of cells (e.g., 10).',...
                'Remove Genes',[1 85],{'0.05'});
            if isempty(answer), return; end
            if iscell(answer)
                a=str2double(answer{1});
                if a>0 && a<intmax
                    fw = gui.gui_waitbar; 
                    sce=sce.selectkeepgenes(1,a);
                    gui.gui_waitbar(fw);
                end
            end
        case 'Remove Genes by Name'        % remove selected genes
            [glist]=gui.i_selectngenes(sce);
            if isempty(glist), return; end
            [y,idx]=ismember(upper(glist),upper(sce.g));
            if ~all(y), error('xxx'); end
            
            %gsorted=sort(sce.g);
            %[idx]=gui.i_selmultidlg(gsorted);
            if isempty(idx), return; end
            if isscalar(idx) && idx==0
                helpdlg('No gene selected.','');
                return;
            end
            %[~,idx]=ismember(sce.g(idx),sce.g);
            
            %{
            answer1 = questdlg(sprintf('Remove %d selected genes?',length(idx)));
            if strcmpi(answer1,'Yes')   
        	    fw = gui.gui_waitbar;                    
                sce.g(idx)=[];
                sce.X(idx,:)=[];
                gui.gui_waitbar(fw);
            else
                return;
            end
            %}
    
            answer1 = questdlg('Remove selected or unselected genes?','', ...
                'Selected', 'Unselected','Selected');
            if isempty(answer1), return; end
            if strcmp(answer1, 'Selected')
        	        fw = gui.gui_waitbar;                    
                    sce.g(idx)=[];
                    sce.X(idx,:)=[];
                    gui.gui_waitbar(fw);
            elseif strcmp(answer1, 'Unselected')
        	        fw = gui.gui_waitbar;                    
                    sce.g=sce.g(idx);
                    sce.X=sce.X(idx,:);
                    gui.gui_waitbar(fw);
            else
                return;
            end

        case 'Remove Mt-genes'        % remove mt-genes
            sce=sce.rmmtgenes;
        case 'Remove Ribosomal Genes'
            sce=sce.rmribosomalgenes;
        case '------------------------------------------------'
            requirerefresh=false;
            return;
        case 'Library Size vs. Mt-reads Ratio'      % mt-ratio vs. library size
            idx=startsWith(sce.g,'mt-','IgnoreCase',true);
            if ~any(idx) 
                disp('No mt genes found.');
                return;
            end
            lbsz=sum(sce.X,1);
            lbsz_mt=sum(sce.X(idx,:),1);
            cj=100*(lbsz_mt./lbsz);
            if issparse(cj), cj=full(cj); end
            ttxtj="mtDNA%";

            ci=sum(sce.X);
            if issparse(ci), ci=full(ci); end
            ttxti="Library Size";
            a=maxk(ci,10);                
            idx=gui.i_setranges3(ci',cj',[0 a(end)],...
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
            idx=gui.i_setranges3(ci',cj',[0 a(end)],...
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
            idx=gui.i_setranges3(ci',cj',[0 a(end)],...
                    [0 b(end)],ttxti,ttxtj);
        case 'QC Metrics in Violin Plots'   % view QC metrics violin
            gui.i_qcviolin(sce.X,sce.g);
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
        if ~isempty(idx) && any(~idx)
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
        else
            requirerefresh=false;
            return;
        end
    end
    guidata(FigureHandle,sce);
end


% function [whitelist]=i_selectwhitelist_DEL(sce)
%     whitelist=[];        
%     answer = questdlg('Genes in whitelist will not be removed. Select whitelist genes?',...
%                 'Whitelist Genes','Yes','No','Cancel','No');
%     switch answer
%         case 'Yes'
%             [gsorted]=gui.i_sortgenenames(sce);
%             if isempty(gsorted), return; end
%             [idx]=gui.i_selmultidlg(gsorted);
%             if isempty(idx), return; end
%             if isscalar(idx) && idx==0, return; end
%             whitelist=gsorted(idx);
%         case 'No'
%             whitelist=[];
%             return;
%         case 'Cancel'
%             whitelist=0;
%             return;
%         otherwise
%             whitelist=0;
%             return;
%     end
% end
