function [requirerefresh] = callback_Select5000Genes(src)


requirerefresh = false;
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if sce.NumGenes<=500
    warndlg('Number of cells is too small.');
    return;
end

spciestag = gui.i_selectspecies(2);
if isempty(spciestag), return; end


prompt = {'Remove Mt-Genes (MT-ND1, MT-ND6, MT-CYB, MT-COI, MT-ATP6, etc)?', ...
    'Remove Ribosomal Genes (RPSA, RPS2, RPS3, RPL3, RPL4, RPLP1, etc)?', ...
    'Remove Genes Without Approved Symbols?', ...
    'Remove Genes Expressed in Less Than m Cells (m = 0.075 or 0.050, 10 or 50)?', ...
    'Keep Top n Highly Variable Genes (HVGs) (n = 5000 or 2000)?'};
dlgtitle = '';
dims = [1, 85];

definput = {'Yes', 'Yes', 'Yes', '0.075', num2str(min([sce.NumGenes,5000]))};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    requirerefresh = false;
    return;
end

fw = gui.gui_waitbar;

if strcmpi(answer{1},'Yes') || strcmpi(answer{1},'Y')
    sce = sce.rmmtgenes;
    disp('Mt-genes removed.');
    requirerefresh = true;
end
if strcmpi(answer{2},'Yes') || strcmpi(answer{2},'Y')
    sce = sce.rmribosomalgenes;
    disp('Ribosomal genes removed.');
    requirerefresh = true;
end

if strcmpi(answer{3},'Yes') || strcmpi(answer{3},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            % T = readtable(fullfile(mfolder, ...
            %     '../resources', 'HGNCBiomart.txt'));
            load(fullfile(mfolder, ...
                '../resources', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '../resources', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx] = ismember(upper(sce.g), upper(ApprovedSymbol));
    a1 = length(sce.g);
    sce.g(~idx) = [];
    sce.X(~idx, :) = [];
    a2 = length(sce.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
    requirerefresh = true;
end

try
    a = str2double(answer{4});
    if a > 0 && a < intmax
        a1 = length(sce.g);
        sce = sce.selectkeepgenes(1, a);
        a2 = length(sce.g);
        fprintf('%d lowly expressed genes found and removed.\n',a1-a2);
        requirerefresh = true;
    end
catch ME
    warning(ME.message);
end

try
    a = str2double(answer{5});
    T = sc_hvg(sce.X, sce.g);
    glist = T.genes(1:min([a, sce.NumGenes]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y)
        errordlg('Runtime error.');
        return;
    end
    sce.g = sce.g(idx);
    sce.X = sce.X(idx, :);
catch ME
    warning(ME.message);
end



guidata(FigureHandle, sce);
gui.gui_waitbar(fw);

end