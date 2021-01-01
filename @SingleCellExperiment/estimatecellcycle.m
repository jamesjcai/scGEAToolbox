function obj = estimatecellcycle(obj,forced)
    if nargin<2, forced=false; end
    if isempty(obj.c_cell_cycle_tx) || forced
        obj.c_cell_cycle_tx=run.SeuratCellCycle(obj.X,obj.g);
        disp('c_cell_cycle_tx added.');
    else
        disp('c_cell_cycle_tx existed.');
    end
end
