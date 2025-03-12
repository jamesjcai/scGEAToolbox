function [parentfig, sce, isui] = gui_getfigsce(src)
   if ~isprop(src, 'Parent') || ~isprop(src.Parent, 'Parent')
       error('Invalid source object: missing parent properties.');
   end

    parentfig = src.Parent.Parent;

    if isempty(parentfig)
        error('Invalid parent figure.');
    end
    if nargout > 1
        sce = guidata(parentfig);
        if ~isa(sce, 'SingleCellExperiment')
            error('sce is not a SingleCellExperiment object.');
        end
    end
    if nargout > 2
        isui = gui.i_isuifig(parentfig);
    end
