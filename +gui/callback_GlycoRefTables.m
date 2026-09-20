function callback_GlycoRefTables(src, ~)
%GUI.CALLBACK_GLYCOREFTABLES  Menu callback: view the curated glyco
%reference tables.
%
%   The +gly analyses are driven by curated maps - which glycan module a
%   lectin reads, which ligand-receptor families need which glycan context,
%   how many glycosylation sites a gene carries. They decide what the
%   analyses can find, so they are worth being able to read directly rather
%   than inferring from results.
%
%   Nothing here touches the data; every entry just loads a table.
%
% See also GLY.LECTINMAP, GLY.LECTINREAGENTS, GLY.NODEANNOT, GLY.GENESETS,
% GLY.ENZONTO.

[FigureHandle, ~] = gui.gui_getfigsce(src);

items = { ...
    'Glycan module to lectin map (drives glyco-lectin communication)'
    'Ligand-receptor glyco-modulation map (drives glyco re-weighting)'
    'Lectin reagents, epitopes and transcripts'
    'Per-gene glycosylation site counts (UniProt)'
    'Curated glycobiology gene sets'
    'GlycoEnzOnto glycosylation pathways'};

[indx, tf] = gui.myListdlg(FigureHandle, items, ...
    'Select a reference table to view:', [], false);
if tf ~= 1 || isempty(indx), return; end

fw = gui.myWaitbar(FigureHandle);
try
    switch indx
        case 1
            T = gly.lectinmap();
        case 2
            T = gly.lectinmap("lr");
        case 3
            T = gly.lectinreagents();
        case 4
            T = gly.nodeannot();
        case 5
            [~, ~, ~, T] = gly.genesets();
        case 6
            [~, ~, ~, T] = gly.enzonto();
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, 'Glyco reference tables');
    return;
end
gui.myWaitbar(FigureHandle, fw);

gui.i_viewtable(T, FigureHandle);
gui.i_exporttable(T, true, 'Tglycoref', 'GlycoReferenceTable', ...
    [], [], FigureHandle);
end
