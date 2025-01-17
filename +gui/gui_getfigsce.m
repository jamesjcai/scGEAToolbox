function [FigureHandle, sce]=gui_getfigsce(src)
   if ~isprop(src, 'Parent') || ~isprop(src.Parent, 'Parent')
       error('Invalid source object: missing parent properties.');
   end

    FigureHandle = src.Parent.Parent;

    if isempty(FigureHandle)
        error('Invalid parent figure.');
    end

    sce = guidata(FigureHandle);

    
    if ~isa(sce, 'SingleCellExperiment')
        error('sce is not a SingleCellExperiment object.');
    end
