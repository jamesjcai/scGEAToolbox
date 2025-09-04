function [sce] = i_selectinfogenes(sce, spciestag, parentfig)
if nargin<3, parentfig = []; end
% ---------------------------------
if nargin < 2 || isempty(spciestag)
    spciestag = gui.i_selectspecies(2, false, parentfig);
end
if isempty(spciestag), return; end

prompt = {'Remove Mt-Genes (MT-ND1, MT-ND6, MT-CYB, MT-COI, MT-ATP6, etc.)?', ...
    'Remove Hemoglobin Genes (HBA1, HBB, Hba-a1, etc.)?', ...
    'Remove Ribosomal Genes (RPSA, RPS2, RPS3, RPL3, RPL4, RPLP1, etc.)?', ...
    'Remove Genes Without Approved Symbols?'};
dlgtitle = '';
dims = [1, 80];

definput = {'Yes', 'Yes', 'Yes', 'Yes'};

if gui.i_isuifig(parentfig)
    answer = gui.myInputdlg(prompt, dlgtitle, definput, parentfig);
else
    answer = inputdlg(prompt, dlgtitle, dims, definput);
end

if isempty(answer), return; end

fw = gui.myWaitbar(parentfig);
if strcmpi(answer{1},'Yes') || strcmpi(answer{1},'Y')
    sce = sce.rmmtgenes;
    disp('Mt-genes removed.');    
end

if strcmpi(answer{2},'Yes') || strcmpi(answer{2},'Y')
    sce = sce.rmhemoglobingenes;
    disp('Hemoglobin genes removed.');    
end

if strcmpi(answer{3},'Yes') || strcmpi(answer{3},'Y')
    sce = sce.rmribosomalgenes;
    disp('Ribosomal genes removed.');    
end

if strcmpi(answer{4},'Yes') || strcmpi(answer{4},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            load(fullfile(mfolder, ...
                '..', 'assets', 'Biomart', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '..', 'assets', 'Biomart', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx0] = ismember(upper(sce.g), upper(ApprovedSymbol));
    a1 = length(sce.g);    
    sce.X(~idx0, :) = [];
    sce.g(~idx0) = [];
    a2 = length(sce.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
end
gui.myWaitbar(parentfig, fw);

% ---------------------------------

