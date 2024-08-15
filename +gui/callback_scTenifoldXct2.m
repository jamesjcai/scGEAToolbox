function callback_scTenifoldXct2(src, ~)

if ~gui.gui_showrefinfo('scTenifoldXct [PMID:36787742]'), return; end

extprogname = 'py_scTenifoldXct2';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[~, cL] = grp2idx(sce.c_batch_id);
[j1, j2, ~, ~] = aaa(cL, sce.c_batch_id);
if isempty(j1) || isempty(j2)
    warndlg('All cells have the same BATCH_ID. Two samples are required.','')
    return; 
end
sce1 = sce.selectcells(j1);
sce2 = sce.selectcells(j2);

if sce1.NumCells < 50 || sce2.NumCells < 50
    [answer] = questdlg('One of samples contains too few cells (n < 50). Continue?','');
    if ~strcmp(answer, 'Yes'), return; end
end


[~, cL] = grp2idx(sce.c_cell_type_tx);
[~, ~, celltype1, celltype2] = aaa(cL, sce.c_cell_type_tx);
if isempty(celltype1) || isempty(celltype2) 
    warndlg('All cells are the same type. Two different cell types are required.','')
    return; 
end

celltype1 = string(celltype1);
celltype2 = string(celltype2);

a1 = sprintf('%s -> %s', celltype1, celltype2);
a2 = sprintf('%s -> %s', celltype2, celltype1);

twosided = false;
[answer] = questdlg('Select direction: Source (ligand) -> Target (receptor)', '', 'Both', a1, a2, 'Both');
switch answer
    case 'Both'
        ct1 = celltype1;
        ct2 = celltype2;
        twosided = true;
    case a1
        ct1 = celltype1;
        ct2 = celltype2;
    case a2
        ct1 = celltype2;
        ct2 = celltype1;
    otherwise
        return;
end

[Tcell, iscomplete] = run.py_scTenifoldXct2(sce1, sce2, ct1, ct2, twosided, wkdir);


if twosided && iscell(Tcell)
    [T1] = Tcell{1};
    [T2] = Tcell{2};
    if ~isempty(T1)
        a = sprintf('%s -> %s', celltype1, celltype2);
        T1 = addvars(T1, repelem(a, size(T1, 1), 1), 'Before', 1);
        T1.Properties.VariableNames{'Var1'} = 'direction';
    end
    if ~isempty(T2)
        a = sprintf('%s -> %s', celltype2, celltype1);
        T2 = addvars(T2, repelem(a, size(T2, 1), 1), 'Before', 1);
        T2.Properties.VariableNames{'Var1'} = 'direction';
    end
    T = [T1; T2];
else
    if ~isempty(Tcell)
        T = Tcell; 
        a = sprintf('%s -> %s', celltype1, celltype2);
        T = addvars(T, repelem(a, size(T, 1), 1), 'Before', 1);
        T.Properties.VariableNames{'Var1'} = 'direction';
    end   
end

% ---- export result
if ~iscomplete
    errordlg('Running time error.', '');
end
if ~isempty(T)
    [b, a] = pkg.i_tempfile("sctendifoldxct");
    writetable(T, b);

    [answer] = questdlg(sprintf('Result has been saved in %s', b), ...
        '', 'Export result...', 'Locate result file...', ...
        'Export result...');
    switch answer
        case 'Locate result file...'
            winopen(a);
            pause(2)
            reshowdlg;
        case 'Export result...'
            gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table');
        otherwise
            winopen(a);
    end
else
    helpdlg('No ligand-receptor pairs are identified.', '');
end


function reshowdlg
    [answerx] = questdlg('Export result to other format?', '');
    switch answerx
        case 'Yes'
            gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table');
        otherwise
            return;
    end
end


end


function [i1, i2, cL1, cL2] = aaa(listitems, ci)
i1 = [];
i2 = [];
cL1 = [];
cL2 = [];
n = length(listitems);
if n < 2, return; end
[indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
    'SelectionMode', 'multiple', ...
    'ListString', listitems, ...
    'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
if tfx == 1
    if numel(indxx) ~= 2
        errordlg('Please select 2 groups');
        return;
    end
    cL1 = listitems(indxx(1));
    cL2 = listitems(indxx(2));
    i1 = ci == cL1;
    i2 = ci == cL2;
end



end


