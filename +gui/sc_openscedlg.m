function [sce, filename] = sc_openscedlg(~, ~, parentfig)
    if nargin<3, parentfig = []; end
    sce = [];
    filename = [];
    list = {'SCE Data File(s) (*.mat)...', ...
        'TXT/TSV/CSV File (*.txt)...', ...
        'Seurat/Rds File (*.rds)...', ...
        'AnnData/H5ad File (*.h5ad)...', ...
        'Loom File (*.loom)...', ...
        '----------------------------------', ...
        '10x Genomics H5 File(s) (*.h5)...', ...
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
        'Import SCE Data from Workspace...', ...
        'Load Example Data...'};

    % preftagname ='scimilmodelpath'
    preftagname ='openscedlgindex';
    defaultindx = getpref('scgeatoolbox', preftagname, length(list));


    if gui.i_isuifig(parentfig)
        [indx, tf] = gui.myListdlg(parentfig, list, 'Select a source');
    else
        [indx, tf] = listdlg('ListString', list, ...
            'SelectionMode', 'single', ...
            'PromptString', {'Select a source:'}, ...
            'ListSize', [220, 310], ...
            'Name', 'Import Data', ...
            'InitialValue', defaultindx);
    end
    if tf ~= 1, return; end
    ButtonName = list{indx};
    setpref('scgeatoolbox', preftagname, indx);
    %         ButtonName = gui.myQuestdlg(parentfig, ('Select Input Data Type', ...
    %                               'SC_SCATTER', ...
    %                               'SCE Data .mat', ...
    %                               '10x Genomics .mtx', ...
    %                               'TSV/CSV .txt', 'SCE Data .mat');
    switch ButtonName
        case 'Simulate Data [PMID:27122128]...'
            try
                [sce] = in_simulatedata;
            catch ME
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
        case 'SCE Data File(s) (*.mat)...'
            %promotesave = false;
            [filenm, pathname] = uigetfile( ...
                {'*.mat', 'SCE Data Files (*.mat)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick SCE Data File(s)','MultiSelect','on');
            % if ~(fname), return; end
            figure(parentfig);
            if isequal(filenm, 0), return; end
            if ~iscell(filenm)
                scefile = fullfile(pathname, filenm);
                try
                    fw = gui.myWaitbar(parentfig);
                    load(scefile, 'sce');
                catch ME
                    gui.myWaitbar(parentfig, fw, true);
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                gui.myWaitbar(parentfig, fw);
            else
                if ~in_multifilesgo, return; end
                answer = gui.myQuestdlg(parentfig, 'Which set operation method to merge data?', ...
                'Merging method', ...
                    {'Intersect', 'Union'}, 'Intersect');
                if ~ismember(answer, {'Union', 'Intersect'}), return; end
                methodtag = lower(answer);
                fw = gui.gui_waitbar_adv;
                try
                    insce = cell(1, length(filenm));
                    for k = 1:length(filenm)
                        filename = fullfile(pathname, filenm{k});
                        if exist(filename,'file')
                            gui.gui_waitbar_adv(fw,k./length(filenm), ...
                                sprintf('Loading %s...', filenm{k}));
                            load(filename, 'sce');
                            insce{k} = sce;
                            metainfo = sprintf("Source: %s", filename);
                            insce{k} = insce{k}.appendmetainfo(metainfo);
                        end
                    end
                        sce = sc_mergesces(insce, methodtag);
                catch ME
                    gui.gui_waitbar_adv(fw);
                    disp(ME.message);
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                gui.gui_waitbar_adv(fw);
            end
        case '10x Genomics MTX File (*.mtx)...'
            %'Matrix/MTX File (*.mtx)...'
            try
                [sce] = gui.i_readmtx;
            catch ME
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
        case 'TXT/TSV/CSV File (*.txt)...'
            [fname, pathname] = uigetfile( ...
                {'*.tsv;*.csv;*.txt', ...
                'TSV/CSV Format Files (*.tsc, *.csv, *.txt)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a tsv/csv/txt format file');
            figure(parentfig);
            if isequal(fname, 0), return; end
            filename = fullfile(pathname, fname);
            fw = gui.myWaitbar(parentfig);
            try
                [X, g] = sc_readtsvfile(filename);
                sce = SingleCellExperiment(X, g);
                metainfo = sprintf("Source: %s", filename);
                sce = sce.appendmetainfo(metainfo);
            catch ME
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
            if isempty(sce)
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, 'File Import Failure.');
                return;
            else
                gui.myWaitbar(parentfig, fw);
            end
        case 'Seurat/Rds File (*.rds)...'
            [fname, pathname] = uigetfile( ...
                {'*.rds', 'Seurat/Rds Format Files (*.rds)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a rds format file');
            figure(parentfig);
            if isequal(fname, 0), return; end
            filename = fullfile(pathname, fname);
            fw = gui.myWaitbar(parentfig);
            try
                [sce] = sc_readrdsfile(filename);
            catch ME
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
            if isempty(sce)
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, 'File Import Failure.');
                return;
            else
                gui.myWaitbar(parentfig, fw);
            end
        case 'AnnData/H5ad File (*.h5ad)...'
            [filenm, pathname] = uigetfile( ...
                {'*.h5ad', 'H5AD Files (*.h5ad)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a H5AD file');
            figure(parentfig);
            if isequal(filenm, 0), return; end
            filename = fullfile(pathname, filenm);
            fw = gui.myWaitbar(parentfig);
            try
                [X, g, b, filename] = sc_readh5adfile(filename);
                if ~isempty(X)
                    sce = SingleCellExperiment(X, g);
                    metainfo = sprintf("Source: %s", filename);
                    sce = sce.appendmetainfo(metainfo);
                    if ~isempty(b), sce.c_cell_id = b; end
                    gui.myWaitbar(parentfig, fw);
                else
                    gui.myWaitbar(parentfig, fw, true);
                    gui.myErrordlg(parentfig, 'File Import Failure.');
                    return;
                end
            catch ME
                disp(ME.message);
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end            
        case '10x Genomics H5 File(s) (*.h5)...'
            [filenm, pathname] = uigetfile( ...
                {'*.h5;*.hdf5', 'HDF5 Files (*.h5)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick 10x Genomics H5 file(s)','MultiSelect','on');
            if isequal(filenm, 0), return; end
            figure(parentfig);
            if iscell(filenm)
                if ~in_multifilesgo, return; end
                answer = gui.myQuestdlg(parentfig, 'Which set operation method to merge data?', ...
                    'Merging method', ...
                    {'Intersect', 'Union'}, 'Intersect');
                if ~ismember(answer, {'Union', 'Intersect'}), return; end
                methodtag = lower(answer);
                fw = gui.gui_waitbar_adv;
                try
                    insce = cell(1, length(filenm));
                    for k = 1:length(filenm)
                        filename = fullfile(pathname, filenm{k});
                        if exist(filename,'file')
                            gui.gui_waitbar_adv(fw,k./length(filenm), ...
                                sprintf('Reading %s...', filenm{k}));
                            [X, g, b, c] = sc_read10xh5file(filename);
                            if ~isempty(X)
                                insce{k} = SingleCellExperiment(X, g);
                                metainfo = sprintf("Source: %s", filename);
                                insce{k} = insce{k}.appendmetainfo(metainfo);                                        
                                if ~isempty(b), insce{k}.c_cell_id = b; end
                                if ~isempty(c)
                                    insce{k}.c_batch_id = c;
                                else
                                    insce{k}.c_batch_id = string(repmat(matlab.lang.makeValidName(filenm{k}), ...
                                        [insce{k}.NumCells, 1]));
                                end
                            end
                        end
                    end                            
                    sce = sc_mergesces(insce, methodtag);
                catch ME
                    gui.gui_waitbar_adv(fw);                    
                    gui.myErrordlg(parentfig, ME.message, ME.identifier);
                    return;
                end
                gui.gui_waitbar_adv(fw);
            else
                filename = fullfile(pathname, filenm);
                if exist(filename,'file')
                    try
                        [X, g, b, c] = sc_read10xh5file(filename);
                        if ~isempty(X)
                            sce = SingleCellExperiment(X, g);
                            metainfo = sprintf("Source: %s", filenm);
                            sce = sce.appendmetainfo(metainfo);
                            if ~isempty(b), sce.c_cell_id = b; end
                            if ~isempty(c), sce.c_batch_id = c; end
                        else
                            return;
                        end
                    catch ME
                        gui.myErrordlg(parentfig, ME.message, ME.identifier);
                        return;
                    end
                end
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
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
        case '10x Genomics ''outs'' Folder...'
            disp('Open 10x Genomics outs folder...');
            selpath = uigetdir;
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
                figure(parentfig);
            end
            if selpath == 0, return; end
            try
                fw = gui.myWaitbar(parentfig);
                [X, g, celllist, ftdone] = sc_read10xdir(selpath);
                gui.myWaitbar(parentfig, fw);
            catch ME
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
            if ~ftdone
                gui.myErrordlg(parentfig, 'Input Error');
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
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
                figure(parentfig);
            end            
            if selpath == 0, return; end
            try
                fw = gui.myWaitbar(parentfig);
                [X, g, celllist, ftdone] = sc_readparsebio(selpath);
                gui.myWaitbar(parentfig, fw);
            catch ME
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
            if ~ftdone
                gui.myErrordlg(parentfig, 'Input Error.');
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
            acc = inputdlg({'Input Number(s) (e.g., GSM3308549-52):'}, ...
                'GEO Accession', [1, 50], {'GSM7855468'});
            if isempty(acc), return; end
            %acc = strtrim(deblank(acc{1}));
            %acc = strrep(acc,' ','');
            acc = regexprep(acc{1},'[^a-zA-Z0-9,;\-]','');
            if isempty(acc) || ~strlength(acc) > 4, return; end
            if contains(acc,'-')
                accx = pkg.i_expandrange(acc);
                if strcmp('Yes', gui.myQuestdlg(parentfig, ...
                        sprintf('Expand accession number series to: %s?', accx)))
                    acc = accx;
                else
                    acc = regexprep(acc,'[^a-zA-Z0-9,;]','');
                end
            end
            if strlength(acc) > 4 && ~isempty(regexp(acc, 'G.+', 'once'))
                accv = unique(strsplit(acc, {',', ';', ' '}), 'stable');
                if length(accv) > 1
                    dmanswer = gui.myQuestdlg(parentfig, 'Download and merge data sets?', ...
                        '', {'Yes', 'Cancel'}, 'Yes');
                    if ~strcmp(dmanswer, 'Yes'), return; end
                    try
                        fw = gui.myWaitbar(parentfig);
                        [sce] = pkg.pipeline_multisamplesmerge(accv, false, parentfig);
                        gui.myWaitbar(parentfig, fw);
                    catch ME
                        gui.myWaitbar(parentfig, fw);
                        gui.myErrordlg(parentfig, ME.message, ME.identifier);
                        return;
                    end
                else                        
                    try
                        fw = gui.myWaitbar(parentfig);
                        [sce] = sc_readgeoaccess(acc);
                        gui.myWaitbar(parentfig, fw);
                    catch ME
                        gui.myWaitbar(parentfig, fw);
                        gui.myErrordlg(parentfig, ME.message, ME.identifier);
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
            definput = {'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4666nnn/GSM4666986/suppl/GSM4666986_BL41_filtered_feature_bc_matrix.h5'};
            answer = inputdlg(prompt, dlgtitle, [1, 100], definput);
            if isempty(answer), return; end
            if ~isempty(answer{1})
                fw = gui.myWaitbar(parentfig);
                files = websave(tempname, answer{1});
                if iscell(files)
                    f = files{1};
                else
                    f = files;
                end
                if isempty(f), error('f1'); end
                fprintf('[X, g, b] = sc_read10xh5file(''%s'');\n', f);
                [X, g, b, c] = sc_read10xh5file(f);
                sce = SingleCellExperiment(X, g);
                metainfo = sprintf("Source: %s", answer{1});
                sce = sce.appendmetainfo(metainfo);
                if ~isempty(b), sce.c_cell_id = b; end
                if ~isempty(c), sce.c_batch_id = c; end
                gui.myWaitbar(parentfig, fw);
            end
        case 'Import SCE Data from Workspace...'
            a = evalin('base', 'whos');
            b = struct2cell(a);
            valididx = ismember(b(4, :), 'SingleCellExperiment');
            if ~valididx
                gui.myHelpdlg(parentfig, 'No SCE in Workspace.');
                return;
            end            
            a = a(valididx);
            b = b(1, valididx);
            [b,idx]=natsort(b);
            a = a(idx);

            if gui.i_isuifig(parentfig)
                [indx, tf] = gui.myListdlg(parentfig, b, 'Select SCE variable:');
            else
                [indx, tf] = listdlg('PromptString', {'Select SCE variable:'}, ...
                    'liststring', b, ...
                    'SelectionMode', 'multiple', 'ListSize', [220, 300]);
            end
            if tf == 1
                if isscalar(indx)
                    sce = evalin('base', a(indx).name);
                elseif length(indx) > 1
                    answer = gui.myQuestdlg(parentfig, 'Which set operation method to merge genes?', ...
                        'Merging method', ...
                        {'Intersect', 'Union'}, 'Intersect');
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
                        fw = gui.myWaitbar(parentfig);
                        sce = sc_mergesces(insce, methodtag);
                    catch ME
                        gui.myWaitbar(parentfig, fw, true);
                        disp(ME.message);
                        gui.myErrordlg(parentfig, ME.message, ME.identifier);
                        return;
                    end
                    gui.myWaitbar(parentfig, fw);
                end
            else
                return;
            end
            %promotesave = false;
        case 'Load Example Data...'
            % if gui.i_isuifig(parentfig)
            %     answerstruced = uiconfirm(parentfig, 'Load processed or raw data?', '', ...
            %         'Options', {'Processed', 'Raw', 'Cancel'}, ...
            %         'DefaultOption', 'Processed', ...
            %         'Icon', 'question');
            % else
                answerstruced = gui.myQuestdlg(parentfig, 'Load processed or raw data?', ...
                    '', {'Processed', 'Raw', 'Cancel'}, 'Processed');
            %end
            if ~(strcmp(answerstruced, 'Processed') || strcmp(answerstruced, 'Raw'))
                return;
            end
            %promotesave = false;
            pw1 = fileparts(mfilename('fullpath'));
            %fprintf('Loading SCE Data File example_data/workshop_example.mat...');
            %tic;
            file1 = fullfile(pw1, '..', 'example_data', 'new_example_sce.mat');
            if ~exist(file1, "file")
                gui.myErrordlg(parentfig, "Example data file does not exist.");
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
    if isa(sce, 'SingleCellExperiment')
        if isstring(sce.c_cell_id) || ischar(sce.c_cell_id) || iscellstr(sce.c_cell_id)
            sce.c_cell_id = matlab.lang.makeUniqueStrings(sce.c_cell_id);
        end
    end
end


function [y] = in_multifilesgo
    [answer]=gui.myQuestdlg(parentfig, 'Multiple files selected. After reading each file, data will be merged. Continue?','');
    switch answer
        case 'Yes'
            y=true;
            return;
    end
    y=false;
end


    function [sce] = in_simulatedata
        sce=[];
        definput = {'3000', '5000'};
        prompt = {'Number of genes:', ...
            'Number of cells:'};
        dlgtitle = 'Simulation Settings';
        answer = inputdlg(prompt, dlgtitle, [1, 50], definput);

        if isempty(answer), return; end
        try
            numgenes = str2double(answer{1});
            numcells = str2double(answer{2});
            assert(isfinite(numgenes) & numgenes==floor(numgenes));
            assert(isfinite(numcells) & numcells==floor(numcells));
            assert((numgenes >= 1) && (numgenes <= 30000));
            assert((numcells >= 1) && (numcells <= 30000));
        catch
            gui.myWarndlg(parentfig, 'Invalid parameter values.');
            return;
        end        
        try
            fw = gui.myWaitbar(parentfig);
            [X] = sc_simudata(numgenes, numcells, 'lun');
            [sce] = SingleCellExperiment(X);
            sce.c_batch_id = string([ones(round(sce.NumCells/2),1);... 
                2*ones(sce.NumCells-round(sce.NumCells/2),1)]);
            %[c, cL] = grp2idx(sce.c);
            gui.myWaitbar(parentfig, fw);
            % guidata(parentfig, sce);
            % in_RefreshAll(src, [], false, false);
        catch ME
            gui.myWaitbar(parentfig, fw,true);
            gui.myErrordlg(parentfig, ME.message, ME.identifier);
        end
    end


