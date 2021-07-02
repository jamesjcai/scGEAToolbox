function [sce]=i_mergeSubCellNames(sce)
    sce.c_cell_type_tx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");
end

% https://www.mathworks.com/help/matlab/ref/erase.html
% https://www.mathworks.com/help/matlab/ref/pattern.html