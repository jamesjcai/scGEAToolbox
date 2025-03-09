function [sce] = pipeline_multisamplesmerge(accv, guiwaitbar)
if nargin < 2, guiwaitbar = true; end
% https://www.cell.com/cell/fulltext/S0092-8674(19)31178-X
% DOI:https://doi.org/10.1016/j.cell.2019.10.028
%{'T11-Apobec-7daytreated','T11-Apobec-Nottreated',...
% 'KPB25Luv-7daytreated','KPB25Luv-Nottreated'}
if nargin < 1
    % GSM4042585,GSM4042586,GSM4042587
    accv = {'GSM4042585', 'GSM4042586', 'GSM4042587', 'GSM4042588'};
end

accv = string(accv);
SCEV = cell(length(accv), 1);
if guiwaitbar
    fw = gui.gui_waitbar_adv;
    gui.gui_waitbar_adv(fw, 0.15);
end
for k = 1:length(accv)
    [sce] = sc_readgeoaccess(strtrim(accv(k)));
    if guiwaitbar
        gui.gui_waitbar_adv(fw, 0.15+0.75*(k / length(accv)));
    end
    if length(accv) > 1
        sce = sce.qcfilterwhitelist(500, 0.15, 5, '');
        try
           sce.c_batch_id=repmat(accv(k),size(sce.c_batch_id));
        catch ME
           warning(ME.message);
        end
    end
    SCEV{k} = sce;
    pause(2);
end

%%
if length(accv) > 1
    sce = sc_mergesces(SCEV,[],true);
else
    sce = SCEV{1};
end

if guiwaitbar, gui.gui_waitbar_adv(fw); end

answerstruced = questdlg('Process merged SCE data (tSNE, clustering, and cell type annotation)?', ...
    '', 'Yes', 'Skip', 'Yes');
if strcmp(answerstruced, 'Yes')
    % [ndim] = gui.i_choose2d3d;
    ndim = 3;
    if ~isempty(ndim)        
        FigureHandle=[];
        [speciestag] = gui.i_selectspecies(2, false, FigureHandle);
        if ~isempty(speciestag)
            if isempty(ndim), return; end
            sce = sce.embedcells('tsne3d', true, true, ndim);
            %k = round(sce.NumCells/100);
            %sce = sce.clustercells(k, 'kmeans', true);
            sce = sce.clustercells([], [], true);
            %sce = pkg.e_celltypes2allclust(sce, speciestag, true);
            sce = sce.assigncelltype(speciestag, false);
        end
    end
end
end

    % 
    % sce = sce.qcfilter;
    % sce = sce.embedcells('tsne3d',true);
    % sce = sce.clustercells([], [], true);
    % sce = sce.assigncelltype(speciestag, false);

