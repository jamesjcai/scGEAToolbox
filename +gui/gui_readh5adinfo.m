function [sce] = gui_readh5adinfo(filename, sce)

xobs = [];
xvar = [];

try
xobs=h5info(filename,'/obs');
catch
end
try
xvar=h5info(filename,'/var');
catch
end
if ~isempty(xobs)
    assert(isequal(xobs.Name, '/obs'))
    n = size(xobs.Datasets, 1);
    obsv = strings(n, 1);
    for k = 1:n
        obsv(k) = string(xobs.Datasets(k).Name);
    end

    switch questdlg('Select Cell Type?')
        case 'Yes'
            [indx, tf] = listdlg('PromptString', ...
                {'Select dataset of cell type:'}, ...
                'SelectionMode', 'single', 'ListString', obsv, ...
                'ListSize', [220 300]);
            if tf == 1
                celltypestr = obsv(indx);
                sce.c_cell_type_tx = h5read(filename, sprintf('/obs/%s', celltypestr));
            end
    end

    switch questdlg('Select Batch ID?')
        case 'Yes'
            [indx, tf] = listdlg('PromptString', ...
                {'Select dataset of batch id:'}, ...
                'SelectionMode', 'single', 'ListString', obsv, ...
                'ListSize', [220 300]);
            if tf == 1                
                sce.c_batch_id = h5read(filename, sprintf('/obs/%s', obsv(indx)));
            end
    end
    
end
if ~isempty(xvar)
    assert(isequal(xvar.Name, '/var'))
    n = size(xvar.Datasets, 1);
    varv = strings(n, 1);
    for k = 1:n
        varv(k) = string(xvar.Datasets(k).Name);        
    end
end
