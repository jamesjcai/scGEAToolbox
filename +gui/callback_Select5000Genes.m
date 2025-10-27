function [requirerefresh] = callback_Select5000Genes(src)

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

oldcn = sce.NumCells;
oldgn = sce.NumGenes;

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


if gui.i_isuifig(FigureHandle)
    answer = gui.myInputdlg(prompt, dlgtitle, definput, FigureHandle);
else
    answer = inputdlg(prompt, dlgtitle, dims, definput);    
end
if isempty(answer)
    requirerefresh = false;
    return;
end


if sce.NumCells*sce.NumGenes<4e8
    sceori = copy(sce);
    % disp('Ready for reversible.');
else
    answer = gui.myQuestdlg(FigureHandle, 'You are about to change the SCE data. This cannot be undone.');
    if ~strcmp(answer, 'Yes'), return; end        
    sceori = [];
end

fw = gui.myWaitbar(FigureHandle);

c = 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmmtgenes;
    disp('Mt-genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmhemoglobingenes;
    disp('Hemoglobin genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = contains(sce.g, 'orf') | contains(sce.g, '-AS') | contains(sce.g, '-as');    
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name contains ''orf'' or ''-AS'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = startsWith(sce.g, 'LINC');    
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name starts with ''LINC'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    sce = sce.rmribosomalgenes;
    disp('Ribosomal genes removed.');
    % requirerefresh = true;
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = find(~cellfun(@isempty, regexp(sce.g,"Gm[0-9][0-9][0-9]")));
    sce.X(idx, :) = [];
    sce.g(idx) = [];    
    a2 = length(sce.g);
    fprintf('%d genes with name starts with ''Gm'' are found and removed.\n',a1-a2);
end


c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    a1 = length(sce.g);
    idx = endsWith(sce.g, 'Rik');    
    sce.X(idx, :) = [];
    sce.g(idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes with name ends with ''Rik'' are found and removed.\n',a1-a2);
end

c = c + 1;
if strcmpi(answer{c},'Yes') || strcmpi(answer{c},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            % T = readtable(fullfile(mfolder, ...
            %     '../assets', 'HGNCBiomart.txt'));
            load(fullfile(mfolder, ...
                '..', 'assets', 'Biomart', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '..', 'assets',  'Biomart', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx] = ismember(upper(sce.g), upper(ApprovedSymbol));
    a1 = length(sce.g);    
    sce.X(~idx, :) = [];
    sce.g(~idx) = [];
    a2 = length(sce.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
    % requirerefresh = true;
end

try
    a = str2double(answer{end-1});
    if a > 0 && a < intmax
        a1 = length(sce.g);
        sce = sce.selectkeepgenes(1, a);
        a2 = length(sce.g);
        fprintf('%d lowly expressed genes found and removed.\n',a1-a2);
        % requirerefresh = true;
    end
catch ME
    warning(ME.message);
end

try
    a = str2double(answer{end});
    % T = sc_hvg(sce.X, sce.g);
    T = sc_splinefit(sce.X, sce.g);
    glist = T.genes(1:min([a, sce.NumGenes]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y)
        gui.myErrordlg(FigureHandle, 'Runtime error.');
        return;
    end    
    sce.X = sce.X(idx, :);
    sce.g = sce.g(idx);
 catch ME
     gui.myWaitbar(FigureHandle, fw,true);
     gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
     return;
 end

try
    sce = sce.qcfilterwhitelist(1000, 0.15, 15, 500, []);
catch ME
    gui.myWaitbar(FigureHandle, fw,true);
    gui.myWarndlg(FigureHandle, ME.message, ME.identifier);
    return;
end

gui.myWaitbar(FigureHandle, fw);

newcn = sce.NumCells;
newgn = sce.NumGenes;

    if newgn==0
        if ~isempty(sceori)            
            gui.myHelpdlg(FigureHandle, "All genes are removed. Opertaion is cancelled.");
            % requirerefresh = false;
            sce = sceori;
        else
            requirerefresh = true;
        end
        return;
    end
    if newcn==0
        if ~isempty(sceori)            
            gui.myHelpdlg(FigureHandle, "All cells are removed. Opertaion is cancelled.");
            % requirerefresh = false;
            sce = sceori;
        else
            requirerefresh = true;
        end
        return;
    end
    if oldcn-newcn==0 && oldgn-newgn==0
        gui.myHelpdlg(FigureHandle, "No cells and genes are removed.");
        % requirerefresh = false;
        return;
    end
    if ~isempty(sceori)
        answer = gui.myQuestdlg(FigureHandle, sprintf('%d genes will be removed; %d cells will be removed.\n[%d genes x %d cells] => [%d genes x %d cells]', ...
                oldgn-newgn, oldcn-newcn, oldgn, oldcn, newgn, newcn),'', ...
                {'Accept Changes', 'Cancel Changes'}, 'Accept Changes');
        if ~strcmp(answer, 'Accept Changes')
            % requirerefresh = false;
            sce = sceori;
        else
            requirerefresh = true;
        end
    else
        requirerefresh = true;
    end
    gui.myGuidata(FigureHandle, sce, src);

end