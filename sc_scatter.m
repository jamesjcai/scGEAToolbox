function sc_scatter(X, genelist, s, c)
% SC_SCATTER
%   SC_SCATTER(X,genelist,s,c) displays circles at the locations specified
%   by s, coordinate of cell embedding, which is an n-by-p matrix
%   specifying coordinates for each cell.
%
%   See also SC_SCATTER_SCE.

if usejava('jvm') && ~feature('ShowFigureWindows')
    error('MATLAB is in a text mode. This function requires a GUI-mode.');
end
promotesave=true;
    if nargin < 1
        list={'SCE Data File (*.mat)...',...
              'Matrix/MTX File (*.mtx)...',...
              'TSV/CSV File (*.txt)...',...
              'Seurat/Rds File (*.rds)...',...
              '10x Genomics Folder...',...
              '----------------------------------',...              
              'Link to GEO mtx.gz File...',... 
              'Link to GEO txt.gz File...',... 
              'GEO Accession Number...',...
              '----------------------------------',...
              'Open Workspace Variable SCE...',...
              'Load Example Data...'};

          [indx,tf] = listdlg('ListString',list,...
            'SelectionMode','single',...
            'PromptString',{'Select an input data type:'},...
            'ListSize',[230,200],...
            'Name','SC_SCATTER');
        if tf~=1, return; end
        ButtonName=list{indx};
        
%         ButtonName = questdlg('Select Input Data Type', ...
%                               'SC_SCATTER', ...
%                               'SCE Data .mat', ...
%                               '10x Genomics .mtx', ...
%                               'TSV/CSV .txt', 'SCE Data .mat');
        switch ButtonName
            case 'SCE Data File (*.mat)...'
                promotesave=false;
                [fname, pathname] = uigetfile( ...
                                              {'*.mat', 'MAT-files (*.mat)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a MAT-file');
                if ~(fname)
                    return
                end
                scefile = fullfile(pathname, fname);
                try
                    fw = gui.gui_waitbar;
                    load(scefile, 'sce');
                catch ME
                    gui.gui_waitbar(fw);
                    errordlg(ME.message);
                    return;
                end
                gui.gui_waitbar(fw);

            case 'Matrix/MTX File (*.mtx)...'
                [fname, pathname] = uigetfile( ...
                                              {'*.mtx', 'MTX Format Files (*.mtx)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a mtx format file');
                if ~(fname), return; end
                prefixstr=extractBefore(fname,max([strfind(fname,'matrix'),1]));                
                matrixmtxfile = fullfile(pathname, fname);
                
                
                featurestxtfile = fullfile(pathname, sprintf('%sfeatures.tsv',prefixstr));
                if ~exist(featurestxtfile, 'file')
                    featurestxtfile = fullfile(pathname, sprintf('%sgenes.tsv',prefixstr));
                end
                if ~exist(featurestxtfile, 'file')
                    featurestxtfile = fullfile(pathname, sprintf('%sfeatures.txt',prefixstr));
                end
                if ~exist(featurestxtfile, 'file')
                    featurestxtfile = fullfile(pathname, sprintf('%sgenes.txt',prefixstr));
                end
                if ~exist(featurestxtfile, 'file')
                    answer = questdlg('Pick features.tsv file?');
                    % error('Cannot find features.tsv')
                    switch answer
                        case 'Yes'
                            [fname2, pathname2] = uigetfile( ...
                                                            {'*.tsv', 'TSV Format Files (*.tsv)'
                                                             '*.*',  'All Files (*.*)'}, ...
                                                            'Pick features.tsv file');
                            if ~(fname2)
                                return;
                            end
                            featurestxtfile = fullfile(pathname2, fname2);
                        otherwise
                            helpdlg('Action Cancelled.','');
                            return;
                    end
                else
                    answer = questdlg(sprintf('Use %s?',featurestxtfile),...
                        'Pick features/genes.tsv file');
                    switch answer
                        case 'Yes'
                        case 'No'
                            helpdlg('Action Cancelled.','');
                            return;
                        otherwise
                            helpdlg('Action Cancelled.','');
                            return;
                    end
                end
                [X, genelist] = sc_readmtxfile(matrixmtxfile, featurestxtfile, [], 2);
                sce = SingleCellExperiment(X, genelist);
            case 'TSV/CSV File (*.txt)...'
                [fname, pathname] = uigetfile( ...
                                              {'*.tsv;*.csv;*.txt', 'TSV/CSV Format Files (*.tsc, *.csv, *.txt)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a tsv/csv/txt format file');
                if ~(fname), return; end
                filename = fullfile(pathname, fname);
                [X, genelist] = sc_readtsvfile(filename);
                sce = SingleCellExperiment(X, genelist);
            case 'Seurat/Rds File (*.rds)...'
                [fname, pathname] = uigetfile( ...
                                              {'*.rds', 'Seurat/Rds Format Files (*.rds)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a rds format file');
                if ~(fname), return; end
                filename = fullfile(pathname, fname);
                fw = gui.gui_waitbar;
                [sce] = sc_readrdsfile(filename);
                gui.gui_waitbar(fw);                
            case '10x Genomics Folder...'
                selpath = uigetdir;
                if selpath==0, return; end
                try
                    fw = gui.gui_waitbar;
                    [X,genelist,~,ftdone]=sc_read10xdir2(selpath);
                    gui.gui_waitbar(fw);
                catch ME
                    gui.gui_waitbar(fw);
                    errordlg(ME.message);
                    return;
                end
                if ~ftdone, errordlg('Input Error'); return; end
                sce = SingleCellExperiment(X, genelist);
            case 'GEO Accession Number...'
                acc=inputdlg({'Input number (e.g., GSM3308545):'},...
                    'GEO Accession',[1 40],{'GSM3308545'});
                if isempty(acc), return; end
                acc=deblank(acc{1});
                if strlength(acc)>4 && ~isempty(regexp(acc,'G.+','once'))
                    try                
                        fw=gui.gui_waitbar;                
                        [sce]=sc_readgeoaccession(acc);
                        gui.gui_waitbar(fw);
                    catch ME
                        gui.gui_waitbar(fw);
                        errordlg(ME.message);
                        return;
                    end
                end
            case 'Link to GEO mtx.gz File...'
                [X,genelist,celllist,ftdone]=gui.i_inputgeolinks;
                if isempty(X) || isempty(genelist) || ~ftdone
                    % errordlg('Input Error');
                    return;
                end
                sce = SingleCellExperiment(X, genelist);
                if ~isempty(celllist) && length(celllist)==sce.NumCells
                    sce.c_cell_id=celllist;
                end
            case 'Link to GEO txt.gz File...'
                    prompt = {'Enter link to counts.txt.gz or counts.csv.gz:'};
                    dlgtitle = 'Input Download Links';
                    dims = [1 100];
                    definput = {'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5350nnn/GSM5350808/suppl/GSM5350808_Fibroblast_young_1wk_Saline_counts.csv.gz'};
                    answer = inputdlg(prompt,dlgtitle,dims,definput);
                    if isempty(answer), return; end
                    if ~isempty(answer{1})
                        tmpd=tempdir;
                        if strcmpi(answer{1}(end-2:end),'.gz')
                            files=gunzip(answer{1},tmpd);
                        elseif strcmpi(answer{1}(end-3:end),'.zip')
                            files=unzip(answer{1},tmpd);
                        else
                            errordlg('File format is not supported.');
                            return;
                        end
                        f=files{1};
                        if isempty(f), error('f1'); end
                        [X,genelist]=sc_readtsvfile(f);
                        sce = SingleCellExperiment(X, genelist);
                    end
            case 'Open Workspace Variable SCE...'
                a=evalin('base','whos');
                b=struct2cell(a);
                valididx=ismember(b(4,:),'SingleCellExperiment');
                if isempty(valididx), return; end
                a=a(valididx);
                [indx,tf]=listdlg('PromptString',{'Select SCE variable:'},...
                    'liststring',b(1,valididx),'SelectionMode','single');
                if tf==1
                    sce=evalin('base',a(indx).name);
                else
                    return;
                end
                promotesave=false;
            case 'Load Example Data...'
                promotesave=false;
                pw1=fileparts(mfilename('fullpath'));
                fprintf('Loading SCE Data File example_data/testSce.mat...');
                tic;
                file1=fullfile(pw1,'example_data','testSce.mat');
                load(file1,'sce');
                fprintf('Done.\n');
                toc;
            otherwise
                return;
        end
    else
        if isa(X, 'SingleCellExperiment')
            sc_scatter_sce(X);
            return;
        end
        if nargin < 4 || isempty(c)
            c = ones(size(X, 2), 1);
        end
        if nargin < 3 || isempty(s)
            s = randn(size(X, 2), 3);
        end
        if nargin < 2 || isempty(genelist)
            genelist = string((1:size(X, 1))');
        end
        sce = SingleCellExperiment(X, genelist, s, c);
    end
    if isempty(sce), return; end
    try
        sc_scatter_sce(sce);
    catch ME
        disp(ME.identifier);
        errordlg(ME.message);
    end
    if promotesave
        labels = {'Save SCE to variable named:'}; 
        vars = {'sce'};
        values = {sce};
        export2wsdlg(labels,vars,values,...
            'Save Data to Workspace');
    end
end
