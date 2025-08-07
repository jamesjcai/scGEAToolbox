function callback_scTenifoldXct2(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

numglist = [1 3000 5000];
memmlist = [16 32 64 128];
neededmem = memmlist(sum(sce.NumGenes > numglist));
[yesgohead, prepare_input_only] = gui.i_memorychecked(neededmem);
if ~yesgohead, return; end
    
extprogname = 'py_scTenifoldXct';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end



[~, cL] = findgroups(string(sce.c_batch_id));
[j1, j2, ~, ~] = aaa(cL, string(sce.c_batch_id), FigureHandle);
if isempty(j1) || isempty(j2)
    gui.myWarndlg(FigureHandle, ['All cells have the same BATCH_ID. ' ...
        'Two samples are required.']);
    return; 
end
sce1 = sce.selectcells(j1);
sce2 = sce.selectcells(j2);

if sce1.NumCells < 50 || sce2.NumCells < 50
    if ~strcmp(gui.myQuestdlg(FigureHandle, 'One of samples contains too few cells (n < 50). Continue?'), 'Yes'), return; end
end


[~, cL] = findgroups(string(sce.c_cell_type_tx));
[~, ~, celltype1, celltype2] = aaa(cL, string(sce.c_cell_type_tx), FigureHandle);
if isempty(celltype1) || isempty(celltype2) 
    gui.myWarndlg(FigureHandle, ['All cells are the same type. ' ...
        'Two different cell types are required.']);
    return; 
end

celltype1 = string(celltype1);
celltype2 = string(celltype2);

a1 = sprintf('%s -> %s', celltype1, celltype2);
a2 = sprintf('%s -> %s', celltype2, celltype1);

twosided = false;
[answer] = gui.myQuestdlg(FigureHandle, 'Select direction: Source (ligand) -> Target (receptor)', ...
    '', {'Both', a1, a2}, 'Both');
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

if ~prepare_input_only
    if ~gui.i_setpyenv([], [], FigureHandle), return; end
end

[Tcell, iscomplete] = run.py_scTenifoldXct2(sce1, sce2, ct1, ct2, twosided, ...
    wkdir, true, prepare_input_only, FigureHandle);

T = [];
if twosided && iscell(Tcell)
    [T1] = Tcell{1};
    [T2] = Tcell{2};
    if istable(T1)
        a = sprintf("%s -> %s", celltype1, celltype2);
        T1 = addvars(T1, repelem(a, height(T1), 1), 'Before', 1);
        T1.Properties.VariableNames{'Var1'} = 'direction';
    end
    if istable(T2)
        a = sprintf("%s -> %s", celltype2, celltype1);
        T2 = addvars(T2, repelem(a, height(T2), 1), 'Before', 1);
        T2.Properties.VariableNames{'Var1'} = 'direction';
    end
    T = [T1; T2];
else
    if ~isempty(Tcell)
        T = Tcell; 
        a = sprintf("%s -> %s", celltype1, celltype2);
        T = addvars(T, repelem(a, height(T), 1), 'Before', 1);
        T.Properties.VariableNames{'Var1'} = 'direction';
    end
end

% ---- export result
if ~prepare_input_only && ~iscomplete
    gui.myErrordlg(FigureHandle, 'Running time error.');
end

if ~isempty(T)

    mfolder = fileparts(mfilename('fullpath'));
    load(fullfile(mfolder, '..', 'assets', 'Ligand_Receptor', ...
         'Ligand_Receptor_more.mat'), 'ligand','receptor');
    % knownpair = false(height(T), 1);
    A = [string(T.ligand) string(T.receptor)];
    B = [ligand receptor];
    [knownpair]= ismember(A, B, 'rows');
    assert(length(knownpair)==height(T));
    T=[T, table(knownpair)];

    % [a, b] = pkg.i_tempdirfile("sctendifoldxct");
    b = matlab.lang.makeValidName(string(datetime));
    b = fullfile(wkdir, b+".txt");
    writetable(T, b);

    T(:,[4 5 6 7 11])=[];
    
    [answer] = gui.myQuestdlg(FigureHandle, sprintf('Result has been saved in %s', b), ...
        '', {'Export result...', 'Locate result file...'}, ...
        'Export result...');
    switch answer
        case 'Locate result file...'
            winopen(a);
            pause(2)
            if strcmp(gui.myQuestdlg(FigureHandle, 'Export result to other format?'), 'Yes')
                gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table',[],[],FigureHandle);
            end
        case 'Export result...'
            gui.i_exporttable(T, false, 'Ttenifldxt2', 'TenifldXt2Table',[],[],FigureHandle);
        otherwise
            winopen(a);
    end
else
    if ~prepare_input_only
        gui.myHelpdlg(FigureHandle, 'No ligand-receptor pairs are identified.');
    else
        if strcmp(gui.myQuestdlg(FigureHandle, 'Input files are prepared successfully. Open working folder?',''), 'Yes')
            winopen(wkdir);
        end
    end    
end

end

function [i1, i2, cL1, cL2] = aaa(listitems, ci, FigureHandle)
    i1 = []; i2 = [];
    cL1 = []; cL2 = [];
    n = length(listitems);
    if n < 2, return; end


   if gui.i_isuifig(FigureHandle)
        [indx, tf] = gui.myListdlg(FigureHandle, listitems, ...
            'Select two groups:', ...
            listitems([n-1, n]));
   else
        [indx, tf] = listdlg('PromptString', {'Select two groups:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
   end
       
    if tf == 1
        if numel(indx) ~= 2
            gui.myErrordlg(FigureHandle, 'Please select 2 groups');            
            return;
        end
        cL1 = listitems(indx(1));
        cL2 = listitems(indx(2));
        i1 = ci == cL1;
        i2 = ci == cL2;
    end
end
