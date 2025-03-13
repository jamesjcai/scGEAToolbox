function [sce] = i_readmtx(parentfig)

if nargin<1, parentfig = []; end

sce = [];

[fname, pathname] = uigetfile( ...
    {'*.mtx', 'MTX Format Files (*.mtx)'; ...
    '*.*', 'All Files (*.*)'}, ...
    'Pick a mtx format file');
figure(parentfig);
if isequal(fname, 0), return; end
prefixstr = extractBefore(fname, max([strfind(fname, 'matrix'), 1]));
matrixmtxfile = fullfile(pathname, fname);


featurestxtfile = fullfile(pathname, sprintf('%sfeatures.tsv', prefixstr));
if ~exist(featurestxtfile, 'file')
    featurestxtfile = fullfile(pathname, sprintf('%sgenes.tsv', prefixstr));
end
if ~exist(featurestxtfile, 'file')
    featurestxtfile = fullfile(pathname, sprintf('%sfeatures.txt', prefixstr));
end
if ~exist(featurestxtfile, 'file')
    featurestxtfile = fullfile(pathname, sprintf('%sgenes.txt', prefixstr));
end
if ~exist(featurestxtfile, 'file')
    answer = gui.myQuestdlg(parentfig, 'Pick features.tsv file?');
    % error('Cannot find features.tsv')
    switch answer
        case 'Yes'
            [fname2, pathname2] = uigetfile( ...
                {'*.tsv', 'TSV Format Files (*.tsv)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick features.tsv file');
            figure(parentfig);
            if ~(fname2)
                gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
                return;
            end
            featurestxtfile = fullfile(pathname2, fname2);
        otherwise
            gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
            return;
    end
else
    answer = gui.myQuestdlg(parentfig, sprintf('Use %s?', featurestxtfile), ...
        'Pick features/genes.tsv file');
    switch answer
        case 'Yes'
        case 'No'
            gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
            return;
        otherwise
            gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
            return;
    end
end


barcodestxtfile = fullfile(pathname, sprintf('%sbarcodes.tsv', prefixstr));
if ~exist(barcodestxtfile, 'file')
    barcodestxtfile = fullfile(pathname, sprintf('%sbarcodes.txt', prefixstr));
end
if ~exist(barcodestxtfile, 'file')
    answer = gui.myQuestdlg(parentfig, 'Pick barcodes.tsv file (optional)?');
    % error('Cannot find features.tsv')
    switch answer
        case 'Yes'
            [fname2, pathname2] = uigetfile( ...
                {'*.tsv', 'TSV Format Files (*.tsv)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick barcodes.tsv file');
            figure(parentfig);
            if ~(fname2)
                barcodestxtfile = [];
            else
                barcodestxtfile = fullfile(pathname2, fname2);
            end
        case 'No'
            barcodestxtfile = [];
        otherwise
            %[X, g] = sc_readmtxfile(matrixmtxfile, featurestxtfile, [], 2);
            %sce = SingleCellExperiment(X, g);
            gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
            return;
    end

else
    answer = gui.myQuestdlg(parentfig, sprintf('Use %s (optional)?', barcodestxtfile), ...
        'Pick barcodes.tsv file');
    switch answer
        case 'Yes'
        case 'No'
            %gui.myHelpdlg(parentfig, ('Action Cancelled.','');
            %return;
            barcodestxtfile = [];
        otherwise
            gui.myHelpdlg(parentfig, 'Action Cancelled.', '');
            return;
    end
end

answer = gui.myQuestdlg(parentfig, sprintf('Matrix file: %s\nFeature file: %s\nBarcode file (optional): %s\nContinne?', ...
    matrixmtxfile, featurestxtfile, barcodestxtfile), 'Confirm File Selection');
if ~strcmp(answer, 'Yes'), return; end
fw = gui.myWaitbar(parentfig);
if exist(matrixmtxfile, 'file') && exist(featurestxtfile, 'file') && exist(barcodestxtfile, 'file')
    [X, g, celllist] = sc_readmtxfile(matrixmtxfile, featurestxtfile, barcodestxtfile, 2);
    sce = SingleCellExperiment(X, g);
    if ~isempty(celllist) && length(celllist) == sce.NumCells
        sce.c_cell_id = celllist;
    end
else
    [X, g] = sc_readmtxfile(matrixmtxfile, featurestxtfile, [], 2);
    sce = SingleCellExperiment(X, g);
end
metainfo = sprintf("Source: %s", matrixmtxfile);
sce = sce.appendmetainfo(metainfo);

gui.myWaitbar(parentfig, fw);
