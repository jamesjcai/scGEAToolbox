function sc_scatter(X, genelist, s, c)
    % SC_SCATTER
    %   SC_SCATTER(X,genelist,s,c) displays circles at the locations specified
    %   by s, coordinate of cell embedding, which is an n-by-p matrix
    %   specifying coordinates for each cell.
    %
    %   See also SC_SCATTER_SCE.

    if nargin < 1
        ButtonName = questdlg('Select Input Data Type', ...
                              'SC_SCATTER', ...
                              'SCE Data .mat', ...
                              '10x Genomics .mtx', ...
                              'TSV/CSV .txt', 'SCE Data .mat');
        switch ButtonName
            case 'SCE Data .mat'
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
                    % rethrow(ME);
                    gui.gui_waitbar(fw);
                    errordlg(ME.message);
                    % rethrow(ME)
                    return
                end
                gui.gui_waitbar(fw);
                sc_scatter_sce(sce);
                return
            case '10x Genomics .mtx'
                [fname, pathname] = uigetfile( ...
                                              {'*.mtx', 'MTX Format Files (*.mtx)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a mtx format file');
                if ~(fname)
                    return
                end
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
                                return
                            end
                            featurestxtfile = fullfile(pathname2, fname2);
                        otherwise
                            return
                    end
                else
                    answer = questdlg(sprintf('Use %s?',featurestxtfile),...
                        'Pick features/genes.tsv file');
                    switch answer
                        case 'Yes'
                        case 'No'
                            return;
                        otherwise
                            return;
                    end
                end
                [X, genelist] = sc_readmtxfile(matrixmtxfile, featurestxtfile, [], 2);
            case 'TSV/CSV txt'
                [fname, pathname] = uigetfile( ...
                                              {'*.tsv;*.csv;*.txt', 'TSV/CSV Format Files (*.tsc, *.csv, *.txt)'
                                               '*.*',  'All Files (*.*)'}, ...
                                              'Pick a tsv/csv/txt format file');
                if ~(fname)
                    return
                end
                filename = fullfile(pathname, fname);
                [X, genelist] = sc_readtsvfile(filename);
            otherwise
                return
        end
%         answer = questdlg('Perform basic QC?');
%         switch answer
%             case 'Yes'
%                 [X, genelist] = sc_qcfilter(X, genelist);
%             case 'No'
%             case 'Cancel'
%                 return
%         end
        sc_scatter(X, genelist);
    else
        if isa(X, 'SingleCellExperiment')
            sc_scatter_sce(X);
            return
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
        try
            sce = SingleCellExperiment(X, genelist, s, c);
        catch
            return
        end
        sc_scatter_sce(sce);
    end
end
