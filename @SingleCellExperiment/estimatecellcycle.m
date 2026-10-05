function obj = estimatecellcycle(obj, forced, methodid)
% methodid: 1 = module scores (sc_cellcyclescore, default)
%           2 = Seurat CellCycleScoring, through R
%           3 = tricycle, through R; also stores the cycle position in
%               radians as the cell attribute "tricycle_position"
if nargin < 3, methodid = 1; end
if nargin < 2, forced = false; end
if isempty(obj.c_cell_cycle_tx) || forced || (isscalar(unique(obj.c_cell_cycle_tx)) && unique(obj.c_cell_cycle_tx)=="undetermined" )
    switch methodid
        case 1
            obj.c_cell_cycle_tx = sc_cellcyclescore(obj.X, obj.g);
        case 2
            obj.c_cell_cycle_tx = run.r_SeuratCellCycle(obj.X, obj.g);
        case 3
            [obj.c_cell_cycle_tx, position] = run.r_tricycle(obj.X, obj.g);
            obj.setCellAttribute("tricycle_position", position);
        otherwise
            error("Unknown methodid %d. Use 1 (module scores), 2 (Seurat), or 3 (tricycle).", methodid);
    end
    disp('SCE.C_CELL_CYCLE_TX added.');
else
    disp('SCE.C_CELL_CYCLE_TX is existing.');
end
end
