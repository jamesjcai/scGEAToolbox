function [requirerefresh, highlightindex] = callback_SelectCellsByQC(src)

mfolder = fileparts(mfilename('fullpath'));
highlightindex = [];
needremove = false;
requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);


% 'SC_QCFILTER (QC Preserves Lowly-expressed Cells/Genes)',...

oldcn = sce.NumCells;
oldgn = sce.NumGenes;

listitems = {'SC_QCFILTER (Basic QC for Cells/Genes)', ...
        'SC_QCFILTER (Enabling Whitelist Genes)', ...
        '------------------------------------------------', ...
        'Remove Empty Genes', ...
        'Remove Genes by Expression', ...
        'Remove Genes by Name or Naming Pattern', ...
        '------------------------------------------------', ...
        '(a) Remove Mt-Genes', ...
        '(b) Remove Hemoglobin Genes', ...
        '(c) Remove Ribosomal Genes', ...
        '(d) Remove Genes Without Approved Symbols', ...
        'Remove Genes (a)+(b)+(c)+(d)', ...
        '------------------------------------------------', ...
        'Remove Cells with No MALAT1 Expression',...
        '------------------------------------------------', ...
        'Library Size vs. Mt-reads Ratio', ...
        'Library Size vs. Number of Genes', ...
        'Abundant lncRNAs vs. Number of Genes'};

       if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, listitems, 'Select Filter');
        else
            [indx, tf] = listdlg('PromptString', {'Select Filter'}, ...
                'SelectionMode', 'single', ...
                'ListString', listitems, ...
                'ListSize', [260, 300]);
        end

    if tf ~= 1
       % requirerefresh = false;
        return;
    end

    switch listitems{indx}

        case {'SC_QCFILTER (Basic QC for Cells/Genes)',...
                'SC_QCFILTER (Enabling Whitelist Genes)'}

            if strcmp(listitems{indx},'SC_QCFILTER (Enabling Whitelist Genes)')
                [whitelist] = gui.i_selectwhitelist(sce, FigureHandle);
                if isnumeric(whitelist)
                    if whitelist==0
                        % requirerefresh=false;
                        return;
                    end
                end
            else
                whitelist = [];
            end


            answer3 = gui.myQuestdlg(FigureHandle, 'Relaxed or Strigent?', ...
                'Cutoff Settings',{'Relaxed (keep more cells/genes)', ...
                'Strigent (remove more cells/genes)'}, ...
                'Strigent (remove more cells/genes)');
            switch answer3
                case 'Relaxed (keep more cells/genes)'
                    definput = {'500', '0.20', '10', '200'};
                case 'Strigent (remove more cells/genes)'
                    definput = {'1000', '0.15', '15', '500'};
                otherwise
                    % requirerefresh = false;
                    return;
            end

            prompt = {'Min reads per cell, i.e., library size (500 or 1000):', ...
                'Max mtDNA ratio per cell (0.15=15% or 0.10=10%):', ...
                'Min nonzero cells per gene (0.01 or 0.05, 10 or 50):', ...
                'Min nonzero genes per cell (200 or 500):'};
            dlgtitle = 'QC Cutoffs';
            dims = [1, 80];
            
            if gui.i_isuifig(FigureHandle)
                answer = gui.myInputdlg(prompt, dlgtitle, definput, FigureHandle);
            else
                answer = inputdlg(prompt, dlgtitle, dims, definput);
            end

            if isempty(answer)
                % requirerefresh = false;
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
                % requirerefresh = false;
                gui.myErrordlg(FigureHandle, 'Invalid input(s).');
                return;
            end

            
            fw = gui.myWaitbar(FigureHandle);

            memerror = false;
            try
                sce = sce.qcfilterwhitelist(libsize, mtratio, ...
                    min_cells_nonzero, numgenes, whitelist);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
                if issparse(sce.X)
                    gui.myWaitbar(FigureHandle, fw, true);
                    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                    % requirerefresh = false;
                    return;
                else
                    memerror = true;
                end
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
                    gui.myWaitbar(FigureHandle, fw, true);
                    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
                    % requirerefresh = false;
                    return;
                end
            end

            if min_cells_nonzero < 1.0
                a = sprintf('%.f%%', min_cells_nonzero*100);
            else
                a = sprintf('%d', min_cells_nonzero);
            end

fprintf('\nCells with more than %.f%% mitochondrial reads or fewer than %d total detected molecules (library size) were excluded. Additionally, cells expressing fewer than %d genes and genes detected in fewer than %s cells were removed from the UMI count matrix, resulting in %d genes across %d cells.\n', ...
            mtratio*100, libsize, numgenes, a, sce.NumGenes, sce.NumCells);



%            fprintf(['\n"We removed cells with more than %.f%% mitochondrial reads and ', ...
%                'with less than %d the total number of detected molecules (library size). '], );
%            fprintf(['\nWe also removed the cells expressing fewer than %d unique genes ', ...
%               'and genes expressed in fewer than %s cells from the non-normalized ', ...
%                'UMI count matrix, resulting in %d genes across %d cells."\n'], ...
%                numgenes, a, sce.NumGenes, sce.NumCells);

            % Based on the visual inspection of the distribution of the detected molecules across the retained cells, we removed cells with fewer than 500 detected transcripts indicating low-quality cells or empty droplets.
            % Genes with few counts (fewer than the 15th percentile based on the distribution of the average gene-wise counts across all cells) were considered uninformative and removed.
            % According to the applied criteria for the quality control of cells and genes, the dataset was finally composed of 12,113 genes and 2,990 cells.

            gui.myWaitbar(FigureHandle, fw);
            %   [Xmajor,Xminor,gmajor,gminor]=pkg.e_makeshadowmat(sce.X,sce.g);
            %   [X1,g1]=pkg.e_shadowmatqc(Xmajor,Xminor,gmajor,gminor);        

        case 'Remove Empty Genes'

            fw = gui.myWaitbar(FigureHandle);
            sce = sce.selectkeepgenes(1, 1);
            gui.myWaitbar(FigureHandle, fw);
            
        case 'Remove Genes by Expression' % remove genes by expression
            answer = inputdlg('Expressed in less than % of cells (e.g., 0.05=5%) or # of cells (e.g., 10).', ...
                'Remove Genes', [1, 50], {'0.05'});
            if isempty(answer), return; end
            if iscell(answer)
                a = str2double(answer{1});
                if a > 0 && a < intmax
                    fw = gui.myWaitbar(FigureHandle);
                    sce = sce.selectkeepgenes(1, a);
                    gui.myWaitbar(FigureHandle, fw);
                end
            end
        case 'Remove Genes by Name or Naming Pattern' % remove selected genes
            answer2 = gui.myQuestdlg(FigureHandle, 'Select genes to be removed.','', ...
                {'Manually Select','Select by Patterns'},'Manually Select');
            switch answer2
                case 'Select by Patterns'
                    prompt = {'Remove Genes With Name Contains ''orf'' or ''-AS'' (C22orf42, C21orf58, etc.)?', ...
                        'Remove Genes With Name Starts With ''LINC'' (LINC01426, LINC01694, etc.)?', ...
                        'Remove Genes With Name Starts With ''Gm'' (Gm12768, Gm13305, etc.)?',...
                        'Remove Genes With Name Ends With ''Rik'' (0610005C13Rik, 0610007C21Ri, etc.)?'};
                    dlgtitle = '';
                    dims = [1, 80];
                    definput = {'Yes', 'Yes', 'Yes', 'Yes'};
                    answer3 = inputdlg(prompt, dlgtitle, dims, definput);
                    if isempty(answer3), return; end
                    c = 1;
                    if strcmpi(answer3{c},'Yes') || strcmpi(answer3{c},'Y')
                        a1 = length(sce.g);
                        idx = contains(sce.g, 'orf') | contains(sce.g, '-AS') | contains(sce.g, '-as');
                        if any(idx)
                            sce.g(idx) = [];
                            sce.X(idx, :) = [];
                        end
                        a2 = length(sce.g);
                        fprintf('%d genes with name contains ''orf'' or ''-AS'' are found and removed.\n',a1-a2);
                    end
                    
                    c = c + 1;
                    if strcmpi(answer3{c},'Yes') || strcmpi(answer3{c},'Y')
                        a1 = length(sce.g);
                        idx = startsWith(sce.g, 'LINC');
                        if any(idx)
                            sce.g(idx) = [];
                            sce.X(idx, :) = [];
                        end
                        a2 = length(sce.g);
                        fprintf('%d genes with name starts with ''LINC'' are found and removed.\n',a1-a2);
                    end

                    c = c + 1;
                    if strcmpi(answer3{c},'Yes') || strcmpi(answer3{c},'Y')
                        a1 = length(sce.g);
                        idx = find(~cellfun(@isempty, regexp(sce.g,"Gm[0-9][0-9][0-9]")));
                        if any(idx)
                            sce.g(idx) = [];
                            sce.X(idx, :) = [];
                        end
                        a2 = length(sce.g);
                        fprintf('%d genes with name starts with ''Gm'' are found and removed.\n',a1-a2);
                    end

                    c = c + 1;
                    if strcmpi(answer3{c},'Yes') || strcmpi(answer3{c},'Y')
                        a1 = length(sce.g);
                        idx = endsWith(sce.g, 'Rik');
                        if any(idx)
                            sce.g(idx) = [];
                            sce.X(idx, :) = [];
                        end
                        a2 = length(sce.g);
                        fprintf('%d genes with name ends with ''Rik'' are found and removed.\n',a1-a2);
                    end

                case 'Manually Select'
                    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
                    if isempty(glist), return; end
                    [y, idx] = ismember(upper(glist), upper(sce.g));
                    if ~all(y), error('Runtime error.'); end
                    if isempty(idx), return; end
                    if isscalar(idx) && idx == 0
                        gui.myHelpdlg(FigureHandle, 'No gene selected.');
                        return;
                    end
                    answer1 = gui.myQuestdlg(FigureHandle, 'Remove selected or unselected genes?', '', ...
                        {'Selected', 'Unselected'}, 'Selected');
                    if isempty(answer1), return; end
                    if strcmp(answer1, 'Selected')
                        fw = gui.myWaitbar(FigureHandle);
                        sce.g(idx) = [];
                        sce.X(idx, :) = [];
                        gui.myWaitbar(FigureHandle, fw);
                    elseif strcmp(answer1, 'Unselected')
                        fw = gui.myWaitbar(FigureHandle);
                        sce.g = sce.g(idx);
                        sce.X = sce.X(idx, :);
                        gui.myWaitbar(FigureHandle, fw);
                    else
                        return;
                    end
                otherwise
            end

        case '(a) Remove Mt-Genes' % remove mt-genes
            sce = sce.rmmtgenes;
        case '(b) Remove Hemoglobin Genes'
            sce = sce.rmhemoglobingenes;            
        case '(c) Remove Ribosomal Genes'
            sce = sce.rmribosomalgenes;
        case '(d) Remove Genes Without Approved Symbols'
            speciestag = gui.i_selectspecies(2, false, FigureHandle);
            if isempty(speciestag)
                % requirerefresh = false;
                return;
            end
            load(fullfile(mfolder, ...
                '..', 'assets', 'Biomart', sprintf('Biomart_%s_genes.mat',speciestag)), 'T');
            ApprovedSymbol = string(T.GeneName);
            [idx] = ~ismember(upper(sce.g), upper(ApprovedSymbol));
            if any(idx)
                answer = gui.myQuestdlg(FigureHandle, sprintf('Remove %d genes lacking approved symbols?', sum(idx)));
                switch answer
                    case 'Yes'
                        fw = gui.myWaitbar(FigureHandle);
                        sce.g(idx) = [];
                        sce.X(idx, :) = [];
                        gui.myWaitbar(FigureHandle, fw);
                    otherwise
                        % requirerefresh = false;
                        return;
                end
            else
                gui.myHelpdlg(FigureHandle, 'No genes found.');
                % requirerefresh = false;
                return;
            end
            % Filter protein-coding genes based on HGNC approval status and remove all non-coding genes and pseudogenes.
            % Filter protein-coding genes with HGNC approved symbols and remove all remaining genes.
            % answer=gui.myQuestdlg(FigureHandle, 'Keep all protein-coding genes with HGNC-approved symbols and remove all remaining genes, such as non-coding genes, pseudogenes, and genes that do not have approved symbols. Continue?','');
        case 'Remove Genes (a)+(b)+(c)+(d)'
            speciestag = gui.i_selectspecies(2, false, FigureHandle);
            if isempty(speciestag)
                % requirerefresh = false;
                return;
            end
            load(fullfile(mfolder, ...
                '..', 'assets',  'Biomart', sprintf('Biomart_%s_genes.mat',speciestag)), 'T');
            ApprovedSymbol = string(T.GeneName);
            sce = sce.rmmtgenes;
            sce = sce.rmhemoglobingenes;
            sce = sce.rmribosomalgenes;
            [idx] = ~ismember(upper(sce.g), upper(ApprovedSymbol));
            if any(idx)
                sce.g(idx) = [];
                sce.X(idx, :) = [];
            end
        case '------------------------------------------------'
            % requirerefresh = false;
            return;
        case 'Remove Cells with No MALAT1 Expression'
            disp('MALAT1 expression indicates cell quality in single-cell RNA sequencing data');
            disp('Zoe A. Clarke, Gary D. Bader');
            disp('bioRxiv 2024.07.14.603469; doi: https://doi.org/10.1101/2024.07.14.603469');
            idx = startsWith(sce.g, 'MALAT1', 'IgnoreCase', true);
            if ~any(idx)
                disp('MALAT1 is not found.');
                return;
            end
            if sum(idx) ~= 1, return; end
            idx = full(sce.X(idx, :) == 0);
            if ~any(idx)
                disp('All cells express MALAT1.');
                needremove = false;
            else
                fprintf('\n%d cells lacking MALAT1 exprssion.\n', sum(idx));
                needremove = true;
                idx = ~idx;
            end
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
            gui.i_qcviolin(sce.X, sce.g, FigureHandle);
            % requirerefresh = false;
            return;
            %         case 11
            %             fw = gui.myWaitbar(FigureHandle);
            %             [Xdecon,contamination]=run.decontX(sce);
            %             sce.X=Xdecon;
            %             guidata(FigureHandle,sce);
            %             gui.myWaitbar(FigureHandle, fw);
            %             figure;
            %             gui.i_stemscatter(sce.s,contamination);
            %             zlabel('Contamination rate')
            %             title('Ambient RNA contamination')
            %             pause(1)
            %             gui.myHelpdlg(FigureHandle, 'Contamination removed.')
            %             requirerefresh=false;
            %             return;
        otherwise
            % requirerefresh = false;
            return;
    end
            

    %if ismember(indx,[7 8 9])
    %if ismember(listitems{idx},{'','',''})
    if needremove
        if issparse(idx), idx = full(idx); end
        if ~isempty(idx) && any(~idx)
            answer = gui.myQuestdlg(FigureHandle, sprintf('Remove or highlight %d cells?', sum(~idx)), ...
                '', {'Remove', 'Highlight', 'Cancel'}, 'Remove');
            switch answer
                case 'Remove'
                    sce = sce.removecells(~idx);
                case 'Highlight'
                    highlightindex = zeros(1, length(idx));
                    highlightindex(~idx) = 1;
                    % requirerefresh = false;
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
        else
            % requirerefresh = false;
            return;
        end
    end

    newcn = sce.NumCells;
    newgn = sce.NumGenes;
    if newgn==0
        gui.myHelpdlg(FigureHandle, "All genes are removed. Opertaion is cancelled.");
        % requirerefresh = false;
        return;
    end
    if newcn==0
        gui.myHelpdlg(FigureHandle, "All cells are removed. Opertaion is cancelled.");
        % requirerefresh = false;
        return;
    end
    if oldcn-newcn==0 && oldgn-newgn==0
        gui.myHelpdlg(FigureHandle, "No cells and genes are removed.");
        % requirerefresh = false;
        return;
    end
    
    answer = gui.myQuestdlg(FigureHandle, sprintf('%d genes will be removed; %d cells will be removed.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
            oldgn-newgn, oldcn-newcn, oldgn, oldcn, newgn, newcn),'', ...
            {'Accept Changes', 'Cancel Changes'}, 'Accept Changes');
    if ~strcmp(answer, 'Accept Changes')
        % requirerefresh = false;
        return;
    end
    requirerefresh = true;
    guidata(FigureHandle, sce);
end


    % function [whitelist]=i_selectwhitelist_DEL(sce)
    %     whitelist=[];
    %     answer = gui.myQuestdlg(FigureHandle, 'Genes in whitelist will not be removed. Select whitelist genes?',...
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
