function [done, CellTypeList, i1, i2, cL1, cL2,...
    outdir] = i_batchmodeprep(sce, prefixtag)

done = false;
CellTypeList=[]; i1=[]; i2=[]; cL1=[]; cL2=[]; outdir=[];

if isscalar(unique(sce.c_cell_type_tx))
    warndlg('Only one cell type or cell type is undetermined.','');
    return;
end

[CellTypeSorted]=pkg.e_sortcatbysize(sce.c_cell_type_tx);
[CellTypeList]=in_selectcelltypes(CellTypeSorted);
if isempty(CellTypeList), return; end

[thisc, clabel]=in_select1class(sce);
if isempty(thisc), return; end
if strcmp(clabel,'Cell Type')
    helpdlg('Cannot select ''Cell Type'' as grouping varialbe.');
    return;
end

%[i1, i2, cL1, cL2] = gui.i_select2grps(sce, false);

[i1, i2, cL1, cL2] = in_twogrpsencoding(thisc);
if isempty(i1) || isempty(i2) || isempty(cL1) || isempty(cL2)
    return;
end
if isscalar(i1) || isscalar(i2), return; end

answer=questdlg('Select a folder to save the outupt Excel files. Continue?','');
if ~strcmp(answer,'Yes'), return; end

outdir = uigetdir;
if ~isfolder(outdir), return; end

needoverwritten=false;
for k=1:length(CellTypeList)
    outfile = sprintf('%s_%s_vs_%s_%s.xlsx', ...
        prefixtag, ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        matlab.lang.makeValidName(string(CellTypeList{k})));
    filesaved = fullfile(outdir, outfile);
    if exist(filesaved,'file')
        needoverwritten=true;
    end
end
if needoverwritten
    answer=questdlg(sprintf('Overwrite existing result file(s) in %s?',outdir),'');    
else
    answer=questdlg(sprintf('Result files will be save in %s. Continue?', outdir), '');
end
if ~strcmp(answer,'Yes'), return; end
done = true;

end


function [CellTypeList]=in_selectcelltypes(CellTypeSorted)
    CellTypeList=[];
    %pause(1);
    %[idx] = gui.i_selmultidlg(CellTypeSorted, CellTypeSorted);
    %if isempty(idx), return; end
    %if idx == 0, return; end
    %CellTypeList = CellTypeSorted(idx);
    [indx2, tf2] = listdlg('PromptString', ...
    {'Select Cell Types:'}, ...
    'SelectionMode', 'multiple', 'ListString', CellTypeSorted, ...
    'InitialValue',1:length(CellTypeSorted), ...
    'ListSize', [220, 300]);
    if tf2 == 1
        CellTypeList = CellTypeSorted(indx2);
    end
end


function [thisc, clabel] = in_select1class(sce, allowunique)
    if nargin < 2, allowunique = false; end
    thisc = [];
    clabel = '';
    
    listitems = {'Current Class (C)'};
    if ~isempty(sce.c_cluster_id)
        if allowunique
            listitems = [listitems, 'Cluster ID'];
        else
            if numel(unique(sce.c_cluster_id)) > 1
                listitems = [listitems, 'Cluster ID'];
            end
        end
    end
    
    if ~isempty(sce.c_cell_cycle_tx)
        if allowunique
            listitems = [listitems, 'Cell Cycle Phase'];
        else
            if numel(unique(sce.c_cell_cycle_tx)) > 1
                listitems = [listitems, 'Cell Cycle Phase'];
            end
        end
    end
    if ~isempty(sce.c_batch_id)
        if allowunique
            listitems = [listitems, 'Batch ID'];
        else
            if numel(unique(sce.c_batch_id)) > 1
                listitems = [listitems, 'Batch ID'];
            end
        end
    end
    
    a = evalin('base', 'whos');
    b = struct2cell(a);
    v = false(length(a), 1);
    for k = 1:length(a)
        if max(a(k).size) == sce.NumCells && min(a(k).size) == 1
            v(k) = true;
        end
    end
    if any(v)
        a = a(v);
        b = b(:, v);
        listitems = [listitems, 'Workspace Variable...'];
    end
    
    [indx2, tf2] = listdlg('PromptString', ...
        {'Select grouping variable:'}, ...
        'SelectionMode', 'single', ...
        'ListString', listitems, ...
        'ListSize', [220, 300]);
    if tf2 == 1
        clabel = listitems{indx2};
        switch clabel
            case 'Current Class (C)'
                thisc = sce.c;
            case 'Cluster ID' % cluster id
                thisc = sce.c_cluster_id;
            case 'Batch ID' % batch id
                thisc = sce.c_batch_id;
            %case 'Cell Type' % cell type
            %    thisc = sce.c_cell_type_tx;
            case 'Cell Cycle Phase' % cell cycle
                thisc = sce.c_cell_cycle_tx;
            case 'Workspace Variable...'
                thisc = i_pickvariable;
        end
    end
    
    function [c] = i_pickvariable
        c = [];
        [indx, tf] = listdlg('PromptString', {'Select variable:'}, ...
            'liststring', b(1, :), ...
            'SelectionMode', 'single', ...
            'ListSize', [220, 300]);
        if tf == 1
            c = evalin('base', a(indx).name);
        end
        %    end
    end
end


function [i1, i2, cL1, cL2]=in_twogrpsencoding(thisc)
    [ci, cLi] = grp2idx(thisc);
    listitems = natsort(string(cLi));
    n = length(listitems);
    if n < 2
        errordlg('Need at least two groups.');
        return;
    end
    [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, ...
        'InitialValue', [n - 1, n], ...
        'ListSize', [220, 300]);
    if tfx == 1
        if numel(indxx) ~= 2
            errordlg('Please select 2 groups');
            return;
        end
        [y1, idx1] = ismember(listitems(indxx(1)), cLi);
        [y2, idx2] = ismember(listitems(indxx(2)), cLi);
        assert(y1 & y2);
        i1 = ci == idx1;
        i2 = ci == idx2;
        cL1 = cLi(idx1);
        cL2 = cLi(idx2);
    else
        i1=[]; i2=[]; cL1=[]; cL2=[];
    end
end
