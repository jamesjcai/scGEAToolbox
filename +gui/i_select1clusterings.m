function [thisc, clabel] = i_select1clusterings(sce)
thisc = [];
clabel = '';

listitems = {''};
methodtagv = fieldnames(sce.struct_cell_clusterings);
for k = 1:length(methodtagv)
    methodtag = methodtagv{k};
    if ~isempty(sce.struct_cell_clusterings.(methodtag))
        listitems = [listitems, methodtag];
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
    listitems = [listitems, 'Customized C...'];
end

listitems(1) = [];
if isempty(listitems)
    helpdlg('No clustering variable is available.', '');
    return;
end

[indx2, tf2] = listdlg('PromptString', ...
    {'Select clustering variable:'}, ...
    'SelectionMode', 'single', 'ListString', listitems);
if tf2 == 1
    clabel = listitems{indx2};
    switch clabel
        case 'Customized C...'
            thisc = i_pickvariable;
        otherwise
            thisc = sce.struct_cell_clusterings.(clabel);
    end
end

    function [c] = i_pickvariable
        c = [];
        [indx, tf] = listdlg('PromptString', {'Select variable:'}, ...
            'liststring', b(1, :), 'SelectionMode', 'single');
        if tf == 1
            c = evalin('base', a(indx).name);
        end
end

end