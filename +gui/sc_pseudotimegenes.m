function sc_pseudotimegenes(sce, t)

    T = sc_hvg(sce.X, sce.g);
    glist = T.genes(1:min([2000, sce.NumGenes]));
    [y, idx] = ismember(glist, sce.g);
    if ~all(y), error('Runtime error.'); end
    sce.g = sce.g(idx);
    sce.X = sce.X(idx, :);
                
    fw = gui.gui_waitbar;
    r = corr(t, sce.X.', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
    gui.gui_waitbar(fw);
    [~, idxp] = maxk(r, 4); % Select top 4 positively correlated genes
    [~, idxn] = mink(r, 3); % Select top 3 negatively correlated genes
    selectedg = sce.g([idxp, idxn]);
    try
        psf1 = figure('WindowStyle', 'modal');
        pkg.i_plot_pseudotimeseries(log2(sce.X+1), ...
            sce.g, t, selectedg);
    catch ME
        if exist('psf1', 'var') && ishandle(psf1)
            close(psf1);
        end
        errordlg(ME.message);
    end

end