function [requirerefresh, highlightindex] = callback_SelectCellsByQC(src)

mfolder = fileparts(mfilename('fullpath'));


needremove = false;
requirerefresh = true;
highlightindex = [];
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
% 'SC_QCFILTER (QC Preserves Lowly-expressed Cells/Genes)',...
listitems = {'SC_QCFILTER (Basic QC for Cells/Genes)', ...
        'Remove Genes by Expression', ...
        'Remove Genes by Name', ...
        'Remove Mt-genes', ...
        'Remove Ribosomal Genes', ...
        '------------------------------------------------', ...
        'Filter Genes with HGNC Approved Symbols', ...
        '------------------------------------------------', ...
        'Library Size vs. Mt-reads Ratio', ...
        'Library Size vs. Number of Genes', ...
        'Abundant lncRNAs vs. Number of Genes', ...
        '------------------------------------------------', ...
        'QC Metrics in Violin Plots'};

    %        '------------- Experimental Options -------------',...
    %        'Remove ambient RNA contamination (R required)'};

    [indx, tf] = listdlg('PromptString', {'Select Filter'}, 'SelectionMode', 'single', ...
        'ListString', listitems, 'ListSize', [250, 300]);
    if tf ~= 1
        requirerefresh = false;
        return;
    end

    switch listitems{indx}

        case 'SC_QCFILTER (Basic QC for Cells/Genes)'            

            answer3 = questdlg('Relaxed or Strigent?', ...
                'Cutoff Settings', 'Relaxed (keep more cells/genes)', ...
                'Strigent (remove more cells/genes)', 'Relaxed (keep more cells/genes)');
            switch answer3
                case 'Relaxed (keep more cells/genes)'
                    definput = {'500', '0.20', '10', '200'};
                case 'Strigent (remove more cells/genes)'
                    definput = {'1000', '0.15', '15', '500'};
                otherwise
                    requirerefresh = false;
                    return;
            end

            prompt = {'Min reads per cell, i.e., library size (500 or 1000):', ...
                'Max mtDNA ratio per cell (0.15=15% or 0.10=10%):', ...
                'Min nonzero cells per gene (0.01 or 0.05, 10 or 50):', ...
                'Min nonzero genes per cell (200 or 500):'};
            dlgtitle = 'QC Cutoffs';
            dims = [1, 65];
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if isempty(answer)
                requirerefresh = false;
                return;
            end
            try
                libsize = str2double(answer{1});
                mtratio = str2double(answer{2});
                min_cells_nonzero = str2double(answer{3});
                numgenes = str2double(answer{4});
                assert((libsize > 0) && (libsize < intmax));
                assert((mtratio >= 0.0) && (mtratio <= 1.0));
                assert((min_cells_nonzero >= 0 && min_cells_nonzero <= 1) || (min_cells_nonzero > 1 && min_cells_nonzero < sce.NumCells));
                assert((numgenes > 0) && (numgenes < intmax));
            catch
                requirerefresh = false;
                errordlg('Invalid input(s).');
                return;
            end

            % [whitelist]=gui.i_selectwhitelist(sce);
            % if isnumeric(whitelist)
            %     if whitelist==0
            %         requirerefresh=false;
            %         return;
            %     end
            % end

            whitelist = [];
            fw = gui.gui_waitbar;

            memerror = false;
            try
                sce = sce.qcfilterwhitelist(libsize, mtratio, ...
                    min_cells_nonzero, numgenes, whitelist);
            catch ME
                % if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))

                if issparse(sce.X)
                    gui.gui_waitbar(fw, true);
                    errordlg(ME.message);
                    requirerefresh = false;
                    return;
                else
                    memerror = true;
                end
            end

            if memerror
                % disp('Making X sparse.');
                if ~isa(sce.X, 'double')
                    [sce.X] = pkg.e_uint2sparse(sce.X);
                else
                    sce.X = sparse(sce.X);
                end

                % disp('Using lite version of QC.');
                % Xobj=refwrap(sce.X);
                % sce.X=[];
                % [g]=pkg.sc_qcfilter_objc(Xobj,sce.g,libsize,mtratio,...
                %         min_cells_nonzero,numgenes);
                % sce.X=Xobj.data;

                %[X,g]=sc_qcfilter_lite(sce.X,sce.g,libsize,mtratio,...
                %        min_cells_nonzero,numgenes);
                %sce=SingleCellExperiment(X,g);
                try
                    sce = sce.qcfilterwhitelist(libsize, mtratio, ...
                        min_cells_nonzero, numgenes, whitelist);
                catch ME
                    gui.gui_waitbar(fw, true);
                    errordlg(ME.message);
                    requirerefresh = false;
                    return;
                end
            end

            if min_cells_nonzero < 1.0
                a = sprintf('%.f%%', min_cells_nonzero*100);
            else
                a = sprintf('%d', min_cells_nonzero);
            end

            fprintf(['"We removed cells with more than %.f%% mitochondrial reads and ', ...
                'with less than %d the total number of detected molecules (library size). '], mtratio*100, libsize);
            fprintf(['We also removed the cells expressing fewer than %d unique genes ', ...
                'and genes expressed in fewer than %s cells from the non-normalized ', ...
                'UMI count matrix, resulting in %d genes across %d cells."\n'], ...
                numgenes, a, sce.NumGenes, sce.NumCells);

            % Based on the visual inspection of the distribution of the detected molecules across the retained cells, we removed cells with fewer than 500 detected transcripts indicating low-quality cells or empty droplets.
            % Genes with few counts (fewer than the 15th percentile based on the distribution of the average gene-wise counts across all cells) were considered uninformative and removed.
            % According to the applied criteria for the quality control of cells and genes, the dataset was finally composed of 12,113 genes and 2,990 cells.

            gui.gui_waitbar(fw);
            %   [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
            %   [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);        

        case 'Remove Genes by Expression' % remove genes by expression
            answer = inputdlg('Expressed in less than % of cells (e.g., 0.05=5%) or # of cells (e.g., 10).', ...
                'Remove Genes', [1, 85], {'0.05'});
            if isempty(answer), return; end
            if iscell(answer)
                a = str2double(answer{1});
                if a > 0 && a < intmax
                    fw = gui.gui_waitbar;
                    sce = sce.selectkeepgenes(1, a);
                    gui.gui_waitbar(fw);
                end
            end
        case 'Remove Genes by Name' % remove selected genes
            [glist] = gui.i_selectngenes(sce);
            if isempty(glist), return; end
            [y, idx] = ismember(upper(glist), upper(sce.g));
            if ~all(y), error('xxx'); end

            %gsorted=sort(sce.g);
            %[idx]=gui.i_selmultidlg(gsorted);
            if isempty(idx), return; end
            if isscalar(idx) && idx == 0
                helpdlg('No gene selected.', '');
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

            answer1 = questdlg('Remove selected or unselected genes?', '', ...
                'Selected', 'Unselected', 'Selected');
            if isempty(answer1), return; end
            if strcmp(answer1, 'Selected')
                fw = gui.gui_waitbar;
                sce.g(idx) = [];
                sce.X(idx, :) = [];
                gui.gui_waitbar(fw);
            elseif strcmp(answer1, 'Unselected')
                fw = gui.gui_waitbar;
                sce.g = sce.g(idx);
                sce.X = sce.X(idx, :);
                gui.gui_waitbar(fw);
            else
                return;
            end        

        case 'Remove Mt-genes' % remove mt-genes
            sce = sce.rmmtgenes;
        case 'Remove Ribosomal Genes'
            sce = sce.rmribosomalgenes;
        case 'Filter Genes with HGNC Approved Symbols'

            % Filter protein-coding genes based on HGNC approval status and remove all non-coding genes and pseudogenes.
            % Filter protein-coding genes with HGNC approved symbols and remove all remaining genes.

            answer=questdlg('Keep all protein-coding genes with HGNC-approved symbols and remove all remaining genes, such as non-coding genes, pseudogenes, and genes that do not have approved symbols. Continue?','');
            switch answer
                case 'Yes'
                    load(fullfile(mfolder, ...
                        '../resources', 'hgnc_coding_genes.mat'), 'ApprovedSymbol');
                    [idx] = ismember(upper(sce.g),upper(ApprovedSymbol));
                    fw = gui.gui_waitbar;
                    sce.g(~idx) = [];
                    sce.X(~idx, :) = [];
                    gui.gui_waitbar(fw);
                otherwise
                    requirerefresh = false;
                    return;
            end
                
        case '------------------------------------------------'
            requirerefresh = false;
            return;
        case 'Library Size vs. Mt-reads Ratio' % mt-ratio vs. library size
            idx = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
            if ~any(idx)
                disp('No mt genes found.');
                return;
            end
            lbsz = sum(sce.X, 1);
            lbsz_mt = sum(sce.X(idx, :), 1);
            cj = 100 * (lbsz_mt ./ lbsz);
            if issparse(cj), cj = full(cj); end
            ttxtj = "mtDNA%";

            ci = sum(sce.X);
            if issparse(ci), ci = full(ci); end
            ttxti = "Library Size";
            a = maxk(ci, 10);
            idx = gui.i_setranges3(ci', cj', [0, a(end)], ...
                [0, 15], ttxti, ttxtj);
            needremove = true;
        case 'Library Size vs. Number of Genes'
            cj = sum(sce.X > 0, 1);
            if issparse(cj), cj = full(cj); end
            ttxtj = "Number of Detected Genes";
            ci = sum(sce.X, 1);
            if issparse(ci), ci = full(ci); end
            ttxti = "Library Size";
            a = maxk(ci, 10);
            b = maxk(cj, 10);
            idx = gui.i_setranges3(ci', cj', [0, a(end)], ...
                [0, b(end)], ttxti, ttxtj);
            needremove = true;
        case 'Abundant lncRNAs vs. Number of Genes' % 'Abundant lncRNAs vs. Number of Genes'
            % remove cells with a high fraction of nuclear lncRNA transcripts
            % (Malat1, Meg3 and Kcnq10t1)
            % https://www.frontiersin.org/articles/10.3389/fncel.2020.00065/full#h3
            cj = sum(sce.X > 0, 1);
            if issparse(cj), cj = full(cj); end
            ttxtj = "Number of Detected Genes";

            idx = matches(upper(sce.g), upper({'Malat1', 'Meg3', 'Kcnq1ot1'}));
            if ~any(idx)
                disp('{Malat1,Meg3,Kcnq1ot1} not found.');
                return;
            end
            ci = sum(sce.X(idx, :), 1);
            if issparse(ci), ci = full(ci); end
            ttxti = "Reads in lncRNAs (Malat1,Meg3,Kcnq1ot1)";

            a = maxk(ci, 10);
            b = maxk(cj, 10);
            idx = gui.i_setranges3(ci', cj', [0, a(end)], ...
                [0, b(end)], ttxti, ttxtj);
            needremove = true;
        case 'QC Metrics in Violin Plots' % view QC metrics violin
            gui.i_qcviolin(sce.X, sce.g);
            requirerefresh = false;
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
            requirerefresh = false;
            return;
    end
            

    %if ismember(indx,[7 8 9])
    %if ismember(listitems{idx},{'','',''})
    if needremove
        if ~isempty(idx) && any(~idx)
            answer = questdlg(sprintf('Remove or highlight %d cells?', sum(~idx)), ...
                '', 'Remove', 'Highlight', 'Cancel', 'Remove');
            switch answer
                case 'Remove'
                    sce = sce.removecells(~idx);
                case 'Highlight'
                    highlightindex = zeros(1, length(idx));
                    highlightindex(~idx) = 1;
                    requirerefresh = false;
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
        else
            requirerefresh = false;
            return;
        end
    end
    guidata(FigureHandle, sce);
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
