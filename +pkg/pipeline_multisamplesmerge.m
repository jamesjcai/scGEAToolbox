function [sce]=pipeline_multisamplesmerge(accv,guiwaitbar)
if nargin<2, guiwaitbar=true; end
% https://www.cell.com/cell/fulltext/S0092-8674(19)31178-X
% DOI:https://doi.org/10.1016/j.cell.2019.10.028
%{'T11-Apobec-7daytreated','T11-Apobec-Nottreated',...
% 'KPB25Luv-7daytreated','KPB25Luv-Nottreated'}
if nargin<1
    % GSM4042585,GSM4042586,GSM4042587
    accv={'GSM4042585','GSM4042586','GSM4042587','GSM4042588'};
end

accv=string(accv);
SCEV=cell(length(accv),1);
if guiwaitbar
    fw=gui.gui_waitbar_adv;
    gui.gui_waitbar_adv(fw,0.15);
end
for k=1:length(accv)
    [sce]=sc_readgeoaccession(accv(k));
    if guiwaitbar
        gui.gui_waitbar_adv(fw,0.15+0.75*(k/length(accv)));
    end
    if length(accv)>1
        sce=sce.qcfilterwhitelist(500,0.15,5,'');
    end
    SCEV{k}=sce;
    pause(2);
end
%%
if length(accv)>1
    sce=sc_mergesces(SCEV);
else
    sce=SCEV{1};
end
if guiwaitbar, gui.gui_waitbar_adv(fw); end

answerstruced=questdlg('Process merged SCE data (tSNE, clustering, and cell type annotation)?',...
    '','Yes','Skip','Yes');
if strcmp(answerstruced,'Yes')
    [speciestag] = gui.i_selectspecies;
    if ~isempty(speciestag)
        [ndim]=gui.i_choose2d3d;
        if isempty(ndim), return; end
        sce = sce.embedcells('tsne', true, true, ndim);
        k=round(sce.NumCells/100);
        sce = sce.clustercells(k, 'kmeans', true);
        sce = pkg.e_celltypes2allclust(sce,speciestag,true);
    end
end
end
