function [requirerefresh, scenew] = callback_Select5000Genes(src)

requirerefresh = false;
scenew = [];

[FigureHandle, sce] = gui.gui_getfigsce(src);

if sce.NumGenes<=500
    gui.myWarndlg(FigureHandle, 'Number of cells is too small.');
    return;
end

spciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(spciestag), return; end


prompt = {'Remove Mt-Genes (MT-ND1, MT-ND6, MT-CYB, MT-COI, MT-ATP6, etc.)?', ...
    'Remove Hemoglobin Genes (HBA1, HBB, Hba-a1, etc.)?', ...
    'Remove Genes With Name Contains ''orf'' or ''-AS'' (C22orf42, C21orf58, etc.)?', ...
    'Remove Genes With Name Starts With ''LINC'' (LINC01426, LINC01694, etc.)?', ...
    'Remove Ribosomal Genes (RPSA, RPS2, RPS3, RPL3, RPL4, RPLP1, etc.)?', ...
    'Remove Genes With Name Starts With ''Gm'' (Gm12768, Gm13305, etc.)?',...
    'Remove Genes With Name Ends With ''Rik'' (0610005C13Rik, 0610007C21Ri, etc.)?',...
    'Remove Genes Without Approved Symbols?', ...
    'Remove Genes Expressed in Less Than m Cells (m = 0.075 or 0.050, 10 or 50)?', ...
    'Keep Top n Highly Variable Genes (HVGs) (n = 5000 or 2000)?'};
dlgtitle = '';
dims = [1, 80];

definput = {'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', '0.075', num2str(min([sce.NumGenes,5000]))};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    requirerefresh = false;
    return;
end

fw = gui.myWaitbar(FigureHandle);
scenew = sce;

c = 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    scenew = scenew.rmmtgenes;
    disp('Mt-genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    scenew = scenew.rmhemoglobingenes;
    disp('Hemoglobin genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(scenew.g);
    idx = contains(scenew.g, 'orf') | contains(scenew.g, '-AS') | contains(scenew.g, '-as');
    scenew.g(idx) = [];
    scenew.X(idx, :) = [];
    a2 = length(scenew.g);
    fprintf('%d genes with name contains ''orf'' or ''-AS'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(scenew.g);
    idx = startsWith(scenew.g, 'LINC');
    scenew.g(idx) = [];
    scenew.X(idx, :) = [];
    a2 = length(scenew.g);
    fprintf('%d genes with name starts with ''LINC'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    scenew = scenew.rmribosomalgenes;
    disp('Ribosomal genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(scenew.g);
    idx = find(~cellfun(@isempty, regexp(scenew.g,"Gm[0-9][0-9][0-9]")));
    scenew.g(idx) = [];
    scenew.X(idx, :) = [];
    a2 = length(scenew.g);
    fprintf('%d genes with name starts with ''Gm'' are found and removed.\n',a1-a2);
end


c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(scenew.g);
    idx = endsWith(scenew.g, 'Rik');
    scenew.g(idx) = [];
    scenew.X(idx, :) = [];
    a2 = length(scenew.g);
    fprintf('%d genes with name ends with ''Rik'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            % T = readtable(fullfile(mfolder, ...
            %     '../resources', 'HGNCBiomart.txt'));
            load(fullfile(mfolder, ...
                '../resources', 'Biomart', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '../resources',  'Biomart', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx] = ismember(upper(scenew.g), upper(ApprovedSymbol));
    a1 = length(scenew.g);
    scenew.g(~idx) = [];
    scenew.X(~idx, :) = [];
    a2 = length(scenew.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
    % requirerefresh = true;
end

try
    a = str2double(answer{end-1});
    if a > 0 && a < intmax
        a1 = length(scenew.g);
        scenew = scenew.selectkeepgenes(1, a);
        a2 = length(scenew.g);
        fprintf('%d lowly expressed genes found and removed.\n',a1-a2);
        % requirerefresh = true;
    end
catch ME
    warning(ME.message);
end

try
    a = str2double(answer{end});
    % T = sc_hvg(scenew.X, scenew.g);
    T = sc_splinefit(scenew.X, scenew.g);
    glist = T.genes(1:min([a, scenew.NumGenes]));
    [y, idx] = ismember(glist, scenew.g);
    if ~all(y)
        gui.myErrordlg(FigureHandle, 'Runtime error.');
        return;
    end
    scenew.g = scenew.g(idx);
    scenew.X = scenew.X(idx, :);
catch ME
    gui.myWaitbar(FigureHandle, fw,true);
    gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
    return;
end

try
    scenew = scenew.qcfilterwhitelist(1000, 0.15, 15, 500, []);
catch ME
    gui.myWaitbar(FigureHandle, fw,true);
    gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
    return;
end

gui.myWaitbar(FigureHandle, fw);
if scenew.NumCells == 0 || scenew.NumGenes==0
    gui.myWarndlg(FigureHandle, 'No cells or genes left. The operation is cancelled');
    requirerefresh = false;
else
    requirerefresh = true;
end

end