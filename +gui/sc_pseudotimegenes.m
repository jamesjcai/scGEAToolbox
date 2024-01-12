function sc_pseudotimegenes(sce, t)


    [K,usehvgs] = gui.i_gethvgnum(sce);
    
    if usehvgs
        % T = sc_hvg(sce.X, sce.g);
        T = sc_splinefit(sce.X, sce.g);
        glist = T.genes(1:min([K, sce.NumGenes]));
        [y, idx] = ismember(glist, sce.g);
        if ~all(y), error('Runtime error.'); end
        sce.g = sce.g(idx);
        sce.X = sce.X(idx, :);
    end

    answer = questdlg('Select method:','','Spearman Correlation','Distance Correlation','Cancel','Spearman Correlation');
    switch answer
        case 'Spearman Correlation'
            fw = gui.gui_waitbar;
            r = corr(t, sce.X.', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
            gui.gui_waitbar(fw);
        case 'Distance Correlation'
            X = log(1+sc_norm(sce.X));
            r = zeros(size(X,1),1);
            n = length(r);
            if n > 100
                fw = gui.gui_waitbar_adv;
                gui.gui_waitbar_adv(fw,100/n);
                for k=1:length(r)
                    if mod(k,100) == 0, gui.gui_waitbar_adv(fw,k/length(r)); end
                    r(k) = pkg.distcorr(X(k,:).', t(:)); % Calculate linear correlation between gene expression profile and T
                end
                gui.gui_waitbar_adv(fw);
            else
                fw = gui.gui_waitbar;
                for k=1:length(r)
                    r(k) = pkg.distcorr(X(k,:).', t(:)); % Calculate linear correlation between gene expression profile and T
                end
                gui.gui_waitbar(fw);
            end
        case 'Cancel'
            return; 
        otherwise
            return;
    end

    [~, idxp] = maxk(r, 10); % Select top 4 positively correlated genes
    % [~, idxn] = mink(r, 3); % Select top 3 negatively correlated genes
    selectedg = sce.g(idxp);
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