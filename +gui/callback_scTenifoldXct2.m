function callback_scTenifoldXct2(src, ~)
gui.gui_showrefinfo('scTenifoldXct [PMID:36787742]');
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
% error('under development.');
[~, cL] = grp2idx(sce.c_batch_id);
[j1, j2, ~, ~] = aaa(cL, sce.c_batch_id);
if isempty(j1) || isempty(j2), return; end


[~, cL] = grp2idx(sce.c_cell_type_tx);
[~, ~, celltype1, celltype2] = aaa(cL, sce.c_cell_type_tx);
if isempty(celltype1) || isempty(celltype2), return; end
sce1 = sce.selectcells(j1);
sce2 = sce.selectcells(j2);

[T, iscomplete] = run.py_scTenifoldXct2(sce1, sce2, celltype1, celltype2);

% ---- export result
if ~iscomplete
    errordlg('Running time error.', '');
end
if ~isempty(T)
    [b, a] = pkg.i_tempfile("sctendifoldxct");
    writetable(T, b);

    answer = questdlg(sprintf('Result has been saved in %s', b), ...
        '', 'Export result...', 'Locate result file...', 'Export result...');
    switch answer
        case 'Locate result file...'
            winopen(a);
            pause(2)
            reshowdlg;
        case 'Export result...'
            gui.i_exporttable(T, false, 'Ttenifldxct', 'TenifldXctTable');

    % "Tcellattrib","CellAttribTable"
    % "Tviolindata","ViolinPlotTable"
    % "Tgenkiglist","GenKIResulTable"
    % "Tmonocleout","MonocleResTable"
    % "Ttenifldxct","TenifldXctTable"

        otherwise
            winopen(a);
    end
else
    helpdlg('No ligand-receptor pairs are identified.', '');
end


% ----- end of export result

end

function [i1, i2, cL1, cL2] = aaa(listitems, ci)
i1 = [];
i2 = [];
cL1 = [];
cL2 = [];
n = length(listitems);
[indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
    'SelectionMode', 'multiple', ...
    'ListString', listitems, ...
    'InitialValue', [n - 1, n]);
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