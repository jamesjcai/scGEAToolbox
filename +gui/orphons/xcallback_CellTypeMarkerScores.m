function [cs, ctselected] = xcallback_CellTypeMarkerScores(src, ~, sce)

answer = questdlg('Calculate a signature score for each cell with respect to a given gene set. Continue?', '');
    if ~strcmp(answer, 'Yes'), return; end

    if nargin < 3
        FigureHandle = src.Parent.Parent;
        sce = guidata(FigureHandle);
    end
    cs = [];
    ctselected = [];
    posg = [];

    actiontype = questdlg('Select a gene set database:', ...
        '', 'PanglaoDB Marker Genes', ...
        'MSigDB Signature Genes', ...
        'DoRothEA TF Target Genes', ...
        'PanglaoDB Marker Genes');

    %[selecteditem] = gui.i_selectcellscore;

    switch actiontype
        case 'PanglaoDB Marker Genes'
            stag = gui.i_selectspecies(2, true);
            if isempty(stag), return; end

            oldpth = pwd;
            pw1 = fileparts(mfilename('fullpath'));
            pth = fullfile(pw1, '..', '+run', 'thirdparty', 'alona_panglaodb');
            cd(pth);

            markerfile = sprintf('marker_%s.mat', stag);
            if exist(markerfile, 'file')
                load(markerfile, 'Tm');
            else
                % Tw=readtable(sprintf('markerweight_%s.txt',stag));
                Tm = readtable(sprintf('markerlist_%s.txt', stag), ...
                    'ReadVariableNames', false, 'Delimiter', '\t');
                % save(markerfile,'Tw','Tm');
            end
            cd(oldpth);

            ctlist = string(Tm.Var1);
            listitems = sort(ctlist);
            [indx, tf] = listdlg('PromptString', ...
                {'Select Class'}, ...
                'SelectionMode', 'single', 'ListString', listitems, 'ListSize', [220, 300]);
            if ~tf == 1, return; end
            ctselected = listitems(indx);
            % idx=find(matches(ctlist,ctselected));
            idx = matches(ctlist, ctselected);
            ctmarkers = Tm.Var2{idx};
            posg = string(strsplit(ctmarkers, ','));
            posg(strlength(posg) == 0) = [];
            ttxt = ctselected;
        case 'MSigDB Signature Genes'
            stag = gui.i_selectspecies(2, true);
            if isempty(stag), return; end

            try
                [posg, ctselected] = gui.i_selectMSigDBGeneSet(stag);
            catch ME
                errordlg(ME.message);
                return;
            end
            ttxt = ctselected;
        case 'DoRothEA TF Target Genes'
            stag = gui.i_selectspecies(2, true);
            if isempty(stag), return; end

            try
                [posg, ctselected] = gui.i_selectTFTargetSet(stag);
            catch ME
                errordlg(ME.message);
                return;
            end
            ttxt = strcat(ctselected, ' activity');
        otherwise
            return;
    end

    if isempty(posg) || isempty(ctselected)
        return;
    end
    % matches(sce.g, posg,'IgnoreCase',true);

    [cs] = gui.e_cellscore(sce, posg);

    disp(['Gene signature scoring - To characterize cells ', ...
        'according to customized/previously reported ', ...
        'gene signature of ___ (Table S1),', ...
        ' gene scores were calculated per cell using ', ...
        'the UCell method [PMID:34285779]/the AddModuleScore function from Seurat.'])

        %     answer = questdlg('Which method?',...
        %     'Select Method', 'AddModuleScore/Seurat', ...
        %     'UCell [PMID: 34285779]','AddModuleScore/Seurat');
        %     switch answer
        %         case 'AddModuleScore/Seurat'
        %             fw=gui.gui_waitbar;
        %             try
        %                 [cs]=sc_cellscore_admdl(sce.X,sce.g,posg);
        %             catch ME
        %                 gui.gui_waitbar(fw,true);
        %                 errordlg(ME.message);
        %             return;
        %             end
        %             gui.gui_waitbar(fw);
        %         case 'UCell [PMID: 34285779]'
        %             %[cs]=run.UCell(sce.X,sce.g,posg);
        %             fw=gui.gui_waitbar;
        %             try
        %                 [cs]=sc_cellscore_ucell(sce.X,sce.g,posg);
        %             catch ME
        %                 gui.gui_waitbar(fw,true);
        %                 errordlg(ME.message);
        %             return;
        %             end
        %             gui.gui_waitbar(fw);
        %         otherwise
        %             return;
        %     end
        %
        %     posg=sort(posg);
        %     fprintf('\n=============\n%s\n-------------\n',ctselected);
        %     for k=1:length(posg)
        %         fprintf('%s\n',posg(k));
        %     end
        %     fprintf('=============\n');
        if nargout == 0
            gui.i_stemscatterfig(sce, cs, posg, ctselected);
            % a=inputdlg('Gene Set Info:','Gene Viewer',[10 50],{char(sce.metadata)});
        end
        %     function i_saveCrossTable(~,~)
        %         gui.i_exporttable(cs,false,'CellScore');
        %     end
        %     function i_geneheatmapx(~,~)
        %         gui.i_geneheatmap(sce,sce.c_cell_type_tx,posg);
        %     end
        answer2 = questdlg(sprintf('Score has been computed.\nCompare the score across cell classes?'), 'Continue?');
        switch answer2
            case 'Yes'

            otherwise
                return;
        end

        [thisc, clabel] = gui.i_select1class(sce);
        if isempty(thisc) % || numel(unique(thisc))==1
            errordlg('Undefined');
            return;
        end

        %     figure('WindowStyle','modal');
        %     pkg.i_violinplot_groupordered(cs,thisc);
        %     ylabel(strrep(ttxt,'_','\_'))
        %     xlabel(clabel);
        gui.i_violinplot(cs, thisc, ttxt);
        xlabel(clabel);

    end
