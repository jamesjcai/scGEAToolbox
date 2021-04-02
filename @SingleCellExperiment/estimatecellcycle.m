function obj = estimatecellcycle(obj,forced,methodid)
    if nargin<3, methodid=1; end
    if nargin<2, forced=false; end
    if isempty(obj.c_cell_cycle_tx) || forced
        switch methodid
            case 1
               obj.c_cell_cycle_tx=sc_cellcyclescoring(obj.X,obj.g);
            case 2
               obj.c_cell_cycle_tx=run.SeuratCellCycle(obj.X,obj.g);
        end
        disp('c_cell_cycle_tx added.');
    else
        disp('c_cell_cycle_tx existed.');
    end
end
