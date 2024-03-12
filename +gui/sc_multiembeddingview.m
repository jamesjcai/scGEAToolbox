function sc_multiembeddingview(sce, embeddingtags, parentfig)
if isempty(embeddingtags)
    embeddingtags=fieldnames(sce.struct_cell_embeddings);
end
    hFig = figure('Visible', 'off');
    hFig.Position(3) = hFig.Position(3) * 1.8;
    for k=1:length(embeddingtags)
        s = sce.struct_cell_embeddings.(embeddingtags{k});
        if size(s,2)>1 && size(s,1)==sce.NumCells
            nexttile
            gui.i_gscatter3(s, sce.c, 1, 1);
            title(embeddingtags{k});
        end
    end
    [px_new] = gui.i_getchildpos(parentfig, hFig);
    if ~isempty(px_new)
        movegui(hFig, px_new);
    else
        movegui(hFig, 'center');
    end
    drawnow;
    hFig.Visible=true;
