function [parentfig, sce, isui] = gui_getfigsce(src)

if nargin < 1 || isempty(src)
    parentfig = [];
    sce = SingleCellExperiment;
    isui = false;
    return;
end

if isa(src, 'matlab.apps.AppBase')
    parentfig = src.UIFigure;
    sce = src.sce;
    isui = true;

elseif isa(src, 'matlab.ui.Figure')
    % A figure is its own parent. Callers that already hold the app window
    % pass it here (gui.i_installtensortoolbox(FigureHandle)); the branch
    % below took src.Parent.Parent, which for a figure is empty, and threw.
    parentfig = src;
    if nargout > 1
        sce = guidata(parentfig);
        if ~isa(sce, 'SingleCellExperiment')
            error('sce is not a SingleCellExperiment object.');
        end
    end
    if nargout > 2
        isui = gui.i_isuifig(parentfig);
    end

else

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

end
