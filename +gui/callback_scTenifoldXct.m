function callback_scTenifoldXct(src,~)
    % import ten.*
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

if ~gui.i_setpyenv
    return; 
end
    
[~,cL]=grp2idx(sce.c_cell_type_tx);
if length(cL)<2, errordlg('Need at least 2 cell types.'); return; end

[indxx,tf2] = listdlg('PromptString',...
    {'Select two cell types:'},...
     'SelectionMode','multiple','ListString',cL);
if tf2==1
    if numel(indxx)~=2
        errordlg('Please select 2 cell types');
        return;
    end
    i1=indxx(1);
    i2=indxx(2);
else
    return;
end
a1=sprintf('Source %s -> Target %s',cL{i1},cL{i2});
a2=sprintf('Source %s -> Target %s',cL{i2},cL{i1});

answer=questdlg('Select direction (ligand->receptor)','',a1,a2,a1);
switch answer
    case a1
        x1=i1; x2=i2;
    case a2
        x1=i2; x2=i1;
    otherwise
        return;
end
idx=sce.c_cell_type_tx==cL{x1} | sce.c_cell_type_tx==cL{x2};
sce=sce.selectcells(idx);

sce.c_batch_id=sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx==cL{x1})="Source";
sce.c_batch_id(sce.c_cell_type_tx==cL{x2})="Target";

species=gui.i_selectspecies(2);
if isempty(species), return; end

try
    fw = gui.gui_waitbar;
    [T]=run.py_scTenifoldXct(sce,species);    
    gui.gui_waitbar(fw);
catch ME
    gui.gui_waitbar(fw,true);
    errordlg(ME.message);
    return;
end
if ~isempty(T)
    gui.i_exporttable(T);
end
end