function [sce] = sc_openscedlg(~, ~)
    sce = [];
    list = {'SCE Data File (*.mat)...', ...
        'TXT/TSV/CSV File (*.txt)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Loom File (*.loom)...', ...
        '----------------------------------', ...
        '10x Genomics H5 File (*.h5)...', ...
        '10x Genomics MTX File (*.mtx)...', ...
        '10x Genomics ''outs'' Folder...', ...
        'Parse Biosciences ''outs'' Folder...', ...
        '----------------------------------', ...
        'Link to GEO h5 File...', ...
        'Link to GEO mtx.gz File...', ...
        'Link to GEO txt.gz File...', ...
        'GEO Accession Number(s)...', ...
        '----------------------------------', ...
        'Simulate Data [PMID:27122128]...',...
        'Import SCE Variable from Workspace...', ...
        'Load Example Data...'};
    [indx, tf] = listdlg('ListString', list, ...
        'SelectionMode', 'single', ...
        'PromptString', {'Select a source:'}, ...
        'ListSize', [230, 295], ...
        'Name', 'Import Data', ...
        'InitialValue', length(list));
    if tf ~= 1, return; end
    ButtonName = list{indx};
    %         ButtonName = questdlg('Select Input Data Type', ...
    %                               'SC_SCATTER', ...
    %                               'SCE Data .mat', ...
    %                               '10x Genomics .mtx', ...
    %                               'TSV/CSV .txt', 'SCE Data .mat');
    switch ButtonName
        case 'Simulate Data [PMID:27122128]...'
            try
                [sce]=in_simulatedata;
            catch ME
                errordlg(ME.message);
                return;
            end
        case 'SCE Data File (*.mat)...'
            promotesave = false;
            [fname, pathname] = uigetfile( ...
                {'*.mat', 'SCE Data Files (*.mat)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a SCE Data File');
            % if ~(fname), return; end
            if isequal(fname, 0), return; end
            scefile = fullfile(pathname, fname);
            try
                fw = gui.gui_waitbar;
                load(scefile, 'sce');
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            gui.gui_waitbar(fw);
        case '10x Genomics MTX File (*.mtx)...'
            %'Matrix/MTX File (*.mtx)...'
            try
                [sce] = gui.i_readmtx;
            catch ME
                errordlg(ME.message);
                return;
            end

        case 'TXT/TSV/CSV File (*.txt)...'
            [fname, pathname] = uigetfile( ...
                {'*.tsv;*.csv;*.txt', 'TSV/CSV Format Files (*.tsc, *.csv, *.txt)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a tsv/csv/txt format file');
            if isequal(fname, 0), return; end
            filename = fullfile(pathname, fname);
            [X, g] = sc_readtsvfile(filename);
            sce = SingleCellExperiment(X, g);
            metainfo = sprintf("Source: %s", filename);
            sce = sce.appendmetainfo(metainfo);
        case 'Seurat/Rds File (*.rds)...'
            [fname, pathname] = uigetfile( ...
                {'*.rds', 'Seurat/Rds Format Files (*.rds)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a rds format file');
            if isequal(fname, 0), return; end
            filename = fullfile(pathname, fname);
            fw = gui.gui_waitbar;
            try
                [sce] = sc_readrdsfile(filename);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            if isempty(sce)
                gui.gui_waitbar(fw, true);
                errordlg('File Import Failure.');
                return;
            else
                gui.gui_waitbar(fw);
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            try
                [X, g, b, filename] = sc_readh5adfile;
                if ~isempty(X)
                    sce = SingleCellExperiment(X, g);
                    metainfo = sprintf("Source: %s", filename);
                    sce = sce.appendmetainfo(metainfo);
                    if ~isempty(b), sce.c_cell_id = b; end
                else
                    return;
                end
            catch ME
                errordlg(ME.message);
                return;
            end
        case '10x Genomics H5 File (*.h5)...'
            try
                [X, g, b, filename] = sc_read10xh5file;
                if ~isempty(X)
                    sce = SingleCellExperiment(X, g);
                    metainfo = sprintf("Source: %s", filename);
                    sce = sce.appendmetainfo(metainfo);
                    if ~isempty(b), sce.c_cell_id = b; end
                else
                    return;
                end
            catch ME
                errordlg(ME.message);
                return;
            end
        case 'Loom File (*.loom)...'
            try
                [X, g, b, filename] = sc_readloomfile;
                if ~isempty(X)
                    sce = SingleCellExperiment(X, g);
                    metainfo = sprintf("Source: %s", filename);
                    sce = sce.appendmetainfo(metainfo);
                    if ~isempty(b), sce.c_cell_id = b; end
                else
                    return;
                end
            catch ME
                errordlg(ME.message);
                return;
            end

        case '10x Genomics ''outs'' Folder...'
            selpath = uigetdir;
            if selpath == 0, return; end
            try
                fw = gui.gui_waitbar;
                [X, g, celllist, ftdone] = sc_read10xdir(selpath);
                gui.gui_waitbar(fw);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            if ~ftdone, errordlg('Input Error');
                return;
            end
            sce = SingleCellExperiment(X, g);
            metainfo = sprintf("Source: %s", selpath);
            sce = sce.appendmetainfo(metainfo);
            if ~isempty(celllist) && length(celllist) == sce.NumCells
                sce.c_cell_id = celllist;
                if isstring(celllist)
                    if all(strlength(celllist) > 17)
                        sce.c_batch_id = extractAfter(sce.c_cell_id, 17);
                        disp('sce.c_batch_id=extractAfter(sce.c_cell_id,17);');
                        disp('Batch IDs assigned.');
                    end
                end
            end
        case 'Parse Biosciences ''outs'' Folder...'
            selpath = uigetdir;
            if selpath == 0, return; end
            try
                fw = gui.gui_waitbar;
                [X, g, celllist, ftdone] = sc_readparsebio(selpath);
                gui.gui_waitbar(fw);
            catch ME
                gui.gui_waitbar(fw, true);
                errordlg(ME.message);
                return;
            end
            if ~ftdone, errordlg('Input Error');
                return;
            end
            sce = SingleCellExperiment(X, g);
            metainfo = sprintf("Source: %s", selpath);
            sce = sce.appendmetainfo(metainfo);
            if ~isempty(celllist) && length(celllist) == sce.NumCells
                sce.c_cell_id = celllist;
                if isstring(celllist)
                    if all(strlength(celllist) > 17)
                        sce.c_batch_id = extractAfter(sce.c_cell_id, 17);
                    end
                end
            end
        case 'GEO Accession Number(s)...'
            acc = inputdlg({'Input Number(s) (e.g., GSM3308547,GSM3308548):'}, ...
                'GEO Accession', [1, 50], {'GSM3308547'});
            if isempty(acc), return; end
            %acc = strtrim(deblank(acc{1}));
            %acc = strrep(acc,' ','');
            acc = regexprep(acc{1},'[^a-zA-Z0-9,;]','');
            if isempty(acc) || ~strlength(acc) > 4, return; end
            if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
                accv = unique(strsplit(acc, {',', ';', ' '}), 'stable');
                if length(accv) > 1
                    dmanswer = questdlg('Download and merge data sets?', ...
                        '', 'Yes', 'Cancel', 'Yes');
                    if ~strcmp(dmanswer, 'Yes'), return; end
                    try
                        fw = gui.gui_waitbar;
                        [sce] = pkg.pipeline_multisamplesmerge(accv, false);
                        gui.gui_waitbar(fw);
                    catch ME
                        gui.gui_waitbar(fw);
                        errordlg(ME.message);
                        return;
                    end
                else                        
                    try
                        fw = gui.gui_waitbar;
                        [sce] = sc_readgeoaccession(acc);
                        gui.gui_waitbar(fw);
                    catch ME
                        gui.gui_waitbar(fw);
                        errordlg(ME.message);
                        return;
                    end
                end
                %                     metainfo=sprintf("Source: %s",acc);
                %                     sce=sce.appendmetainfo(metainfo);
            end
        case {'Link to GEO mtx.gz File...', 'Link to GEO txt.gz File...'}
            if contains(ButtonName, 'mtx')
                [X, g, celllist, ftdone, answer1] = gui.i_inputgeolink_mtx;
            else
                [X, g, celllist, ftdone, answer1] = gui.i_inputgeolink_txt;
            end
            if isempty(X) || isempty(g) || ~ftdone
                return;
            end
            sce = SingleCellExperiment(X, g);
            metainfo = sprintf("Source: %s", answer1);
            sce = sce.appendmetainfo(metainfo);
            if ~isempty(celllist) && length(celllist) == sce.NumCells
                sce.c_cell_id = celllist;
            end
        case 'Link to GEO h5 File...'
            prompt = {'Enter link to .h5 file:'};
            dlgtitle = 'Input Download Links';
            dims = [1, 100];
            definput = {'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4666nnn/GSM4666986/suppl/GSM4666986_BL41_filtered_feature_bc_matrix.h5'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if isempty(answer), return; end
            if ~isempty(answer{1})
                fw = gui.gui_waitbar;
                files = websave(tempname, answer{1});
                if iscell(files)
                    f = files{1};
                else
                    f = files;
                end
                if isempty(f), error('f1'); end
                fprintf('[X,g,b]=sc_read10xh5file(''%s'');\n', f);
                [X, g, b] = sc_read10xh5file(f);
                sce = SingleCellExperiment(X, g);
                metainfo = sprintf("Source: %s", answer{1});
                sce = sce.appendmetainfo(metainfo);
                if ~isempty(b), sce.c_cell_id = b; end
                gui.gui_waitbar(fw);
            end
        case 'Import SCE Variable from Workspace...'
            a = evalin('base', 'whos');
            b = struct2cell(a);
            valididx = ismember(b(4, :), 'SingleCellExperiment');
            if isempty(valididx)
                helpdlg('No SCE in the Workspace.', '');
                return;
            end
            a = a(valididx);
            [indx, tf] = listdlg('PromptString', {'Select SCE variable:'}, ...
                'liststring', b(1, valididx), 'SelectionMode', 'multiple');
            if tf == 1
                if length(indx) == 1
                    sce = evalin('base', a(indx).name);
                elseif length(indx) > 1
                    answer = questdlg('Which set operation method to merge genes?', 'Merging method', ...
                        'Intersect', 'Union', 'Intersect');
                    if ~ismember(answer, {'Union', 'Intersect'}), return; end
                    methodtag = lower(answer);
                    try
                        insce = cell(1, length(indx));
                        s = "";
                        for k = 1:length(indx)
                            insce{k} = evalin('base', a(indx(k)).name);
                            s = sprintf('%s,%s', s, a(indx(k)).name);
                        end
                        s = s(2:end);
                        fprintf('>> sce=sc_mergesces({%s},''%s'');\n', s, methodtag);
                        fw = gui.gui_waitbar;
                        sce = sc_mergesces(insce, methodtag);
                    catch ME
                        gui.gui_waitbar(fw, true);
                        errordlg(ME.message);
                        return;
                    end
                    gui.gui_waitbar(fw);
                end
            else
                return;
            end
            promotesave = false;
        case 'Load Example Data...'
            answerstruced = questdlg('Load processed or raw data?', ...
                '', 'Processed', 'Raw', 'Cancel', 'Processed');
            if ~(strcmp(answerstruced, 'Processed') || strcmp(answerstruced, 'Raw'))
                return;
            end
            promotesave = false;
            pw1 = fileparts(mfilename('fullpath'));
            %fprintf('Loading SCE Data File example_data/workshop_example.mat...');
            %tic;
            file1 = fullfile(pw1, '..', 'example_data', 'workshop_example.mat');
            if ~exist(file1, "file")
                errordlg("Example data file does not exist.");
                return;
            end
            load(file1, 'sce');
            if strcmp(answerstruced, 'Raw')
                orisce = sce;
                sce = SingleCellExperiment(sce.X, sce.g);
                sce.c_batch_id = orisce.c_batch_id;
                sce.c_cell_id = orisce.c_cell_id;
                sce.metadata = orisce.metadata;
                clearvars orisce
            end
            %fprintf('Done.\n');
            %toc;
        case '----------------------------------'
            return;
        otherwise
            return;
    end
end


    function [sce]=in_simulatedata
        sce=[];
        definput = {'3000', '5000'};
        prompt = {'Number of genes:', ...
            'Number of cells:'};
        dlgtitle = 'Simulation Settings';
        dims = [1, 55];
        answer = inputdlg(prompt, dlgtitle, dims, definput);

        if isempty(answer), return; end
        try
            numgenes = str2double(answer{1});
            numcells = str2double(answer{2});
            assert(isfinite(numgenes) & numgenes==floor(numgenes));
            assert(isfinite(numcells) & numcells==floor(numcells));
            assert((numgenes >= 1) && (numgenes <= 30000));
            assert((numcells >= 1) && (numcells <= 30000));
        catch
            errordlg('Invalid parameter value(s).');
            return;
        end
        
            try
                fw = gui.gui_waitbar;
                [X] = sc_simudata(numgenes, numcells,'lun');
                [sce] = SingleCellExperiment(X);
                %[c, cL] = grp2idx(sce.c);
                gui.gui_waitbar(fw);
                % guidata(FigureHandle, sce);
                % in_RefreshAll(src, [], false, false);
            catch ME
                gui.gui_waitbar(fw,true);
                errordlg(ME.message);
            end
    end


