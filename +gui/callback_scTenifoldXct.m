function callback_scTenifoldXct(src, ~)

if ~gui.gui_showrefinfo('scTenifoldXct [PMID:36787742]', FigureHandle), return; end

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


if ~prepare_input_only
    if ~gui.i_setpyenv, return; end
end

[thisc, clabel] = gui.i_select1class(sce, false, 'Select grouping variable (cell type):', 'Cell Type');
if isempty(thisc), return; end

if ~strcmp(clabel, 'Cell Type')
    if ~strcmp(gui.myQuestdlg(FigureHandle, 'You selected grouping varible other than ''Cell Type''. Continue?'), 'Yes'), return; end
end

[c, cL] = grp2idx(thisc);
[idx] = gui.i_selmultidlg(cL, [], FigureHandle);
if isempty(idx), return; end
if numel(idx) < 2
    gui.myWarndlg(FigureHandle, 'Need at least 2 cell groups to perform cell-cell interaction analysis.');
    return;
end
if numel(idx) ~= 2
    gui.myWarndlg(FigureHandle, sprintf('Need only 2 cell groups to perform cell-cell interaction analysis. You selected %d.', ...
        numel(idx)));
    return;
end

i1 = idx(1);
i2 = idx(2);

%{
[~,cL]=grp2idx(sce.c_cell_type_tx);
if length(cL)<2, errordlg('Need at least 2 cell types.'); return; end

[indxx,tf2] = listdlg('PromptString',...
    {'Select two cell types:'},...
    'SelectionMode','multiple','ListString',cL, 'ListSize', [220, 300]);
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
%}


a1 = sprintf('%s -> %s', cL{i1}, cL{i2});
a2 = sprintf('%s -> %s', cL{i2}, cL{i1});

twosided = false;
answer = gui.myQuestdlg(FigureHandle, 'Select direction: Source (ligand) -> Target (receptor)', '', ...
    {'Both', a1, a2}, 'Both');
switch answer
    case 'Both'
        x1 = i1;
        x2 = i2;
        twosided = true;
    case a1
        x1 = i1;
        x2 = i2;
    case a2
        x1 = i2;
        x2 = i1;
    otherwise
        return;
end

%{
idx=sce.c_cell_type_tx==cL{x1} | sce.c_cell_type_tx==cL{x2};
sce=sce.selectcells(idx);

sce.c_batch_id=sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx==cL{x1})="Source";
sce.c_batch_id(sce.c_cell_type_tx==cL{x2})="Target";
%}


sce.c_batch_id = thisc;
sce.c_batch_id(c == x1) = "Source";
sce.c_batch_id(c == x2) = "Target";
sce.c_cell_type_tx = string(cL(c));

% idx=thisc==cL{x1} | thisc==cL{x2};
idx = c == x1 | c == x2;
sce = sce.selectcells(idx);


%sce.c_batch_id(thisc==cL{x1})="Source";
%sce.c_batch_id(thisc==cL{x2})="Target";
T = [];
try
    if twosided
        [Tcell] = run.py_scTenifoldXct(sce, cL{x1}, cL{x2}, true, ...
            wkdir, true, prepare_input_only);
        if ~isempty(Tcell)
            [T1] = Tcell{1};
            [T2] = Tcell{2};
            if ~isempty(T1)
                a = sprintf('%s -> %s', cL{x1}, cL{x2});
                T1 = addvars(T1, repelem(a, height(T1), 1), 'Before', 1);
                T1.Properties.VariableNames{'Var1'} = 'direction';
            end
            if ~isempty(T2)
                a = sprintf('%s -> %s', cL{x2}, cL{x1});
                T2 = addvars(T2, repelem(a, height(T2), 1), 'Before', 1);
                T2.Properties.VariableNames{'Var1'} = 'direction';
            end
            T = [T1; T2];
        end
    else
        [T] = run.py_scTenifoldXct(sce, cL{x1}, cL{x2}, false, wkdir, ...
            true, prepare_input_only);
        %T=readtable('output1.txt');
        if ~isempty(T)
            a = sprintf('%s -> %s', cL{x1}, cL{x2});
            T = addvars(T, repelem(a, height(T), 1), 'Before', 1);
            T.Properties.VariableNames{'Var1'} = 'direction';
        end
    end
catch ME
    errordlg(ME.message);
    return;
end

if ~isempty(T)
    mfolder = fileparts(mfilename('fullpath'));
    load(fullfile(mfolder, '..', 'resources', 'Ligand_Receptor', ...
         'Ligand_Receptor_more.mat'), 'ligand','receptor');
    % knownpair = false(height(T), 1);
    A = [string(T.ligand) string(T.receptor)];
    B = [ligand receptor];
    [knownpair]= ismember(A, B, 'rows');
    assert(length(knownpair)==height(T));

    T=[T, table(knownpair)];
    %T(:,[4 5 6 7 11])=[];
    
    outfile = fullfile(wkdir,"outfile.csv");

    if isfile(outfile)
        answerx = gui.myQuestdlg(FigureHandle, sprintf('Overwrite %s? Select No to save in a temporary file.', outfile));
    else
        answerx = 'Yes';
    end
    if isempty(wkdir) || ~isfolder(wkdir) || ~strcmp(answerx, 'Yes')
        [a, b] = pkg.i_tempdirfile("sctendifoldxct");
        writetable(T, b);
   
        answer = gui.myQuestdlg(FigureHandle, sprintf('Result has been saved in %s', b), ...
            '', {'Export result...', 'Locate result file...'}, 'Export result...');
        switch answer
            case 'Locate result file...'
                winopen(a);
                pause(2)
                if strcmp(gui.myQuestdlg(FigureHandle, 'Export result to other format?'), 'Yes')
                    gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable');
                end
            case 'Export result...'
                gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable');
            otherwise
                winopen(a);
        end
    else
        writetable(T, outfile);
        if strcmp(gui.myQuestdlg(FigureHandle, sprintf('Result has been saved in %s. Open the working folder?', outfile)), 'Yes')
            winopen(wkdir);
        end
    end
else
    if ~prepare_input_only
        gui.myHelpdlg(FigureHandle, 'No ligand-receptor pairs are identified.', '');
    else
        if strcmp(gui.myQuestdlg(FigureHandle, 'Input files are prepared successfully. Open working folder?',''), 'Yes')
            winopen(wkdir);
        end
    end
end

end