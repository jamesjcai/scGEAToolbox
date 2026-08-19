function callback_ViewCellAttributeTable(src, ~)
% CALLBACK_VIEWCELLATTRIBUTETABLE - Show the per-cell attribute table.
%
%   gui.callback_ViewCellAttributeTable(src)
%
% Displays one row per cell with the standard fields (cell id, batch, cluster,
% cell type, cell cycle) followed by everything in sce.list_cell_attributes,
% which is where per-cell metadata imported from Seurat/RDS is stored.
%
% The viewer also exports the table to CSV, Excel, MAT or the workspace, so this
% backs the combined View/Export Cell Attribute Table menu item.
%
% see also: pkg.i_makeattributestable, gui.TableViewerApp, sc_readrdsfile

[parentfig, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells == 0
    gui.myWarndlg(parentfig, 'No data loaded.');
    return;
end

T = pkg.i_makeattributestable(sce);
gui.TableViewerApp(T, parentfig, "CellAttribTable");
end
