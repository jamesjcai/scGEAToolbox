function [requirerefresh] = callback_Select5000Genes(src)

requirerefresh = false;
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if sce.NumGenes<=5000
    return;
end

prompt = {'Remove mt-genes?', ...
    'Remove ribosomal genes?', ...
    'Remove genes expressed less than p cells (e.g., p=0.075)?', ...
    'Select top n HVGs (e.g., n=5000)?'};
dlgtitle = '';
dims = [1, 65];
definput = {'Y', 'Y', '0.075','5000'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer)
    requirerefresh = false;
    return;
end

fw = gui.gui_waitbar;

if strcmpi(answer{1},'Y')
    sce = sce.rmmtgenes; 
    requirerefresh = true;
end
if strcmpi(answer{2},'Y')
    sce = sce.rmribosomalgenes;
    requirerefresh = true;
end

try
    a = str2double(answer{3});
    if a > 0 && a < intmax
        
        sce = sce.selectkeepgenes(1, a);
        requirerefresh = true;

    end
catch ME
    warning(ME.message);
end

try
    a = str2double(answer{4});
    T = sc_hvg(sce.X, sce.g);
    glist = T.genes(1:min([a, sce.NumGenes]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y)
        errordlg('Runtime error.');
        return;
    end
    sce.g = sce.g(idx);
    sce.X = sce.X(idx, :);
catch ME
    warning(ME.message);
end
guidata(FigureHandle, sce);
gui.gui_waitbar(fw);

end