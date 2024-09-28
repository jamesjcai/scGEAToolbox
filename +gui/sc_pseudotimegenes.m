function sc_pseudotimegenes(sce, t, parentfig)

if nargin<3, parentfig=[]; end


    [K,usehvgs] = gui.i_gethvgnum(sce);
    
    if usehvgs
        T = sc_splinefit(sce.X, sce.g);
        glist = T.genes(1:min([K, sce.NumGenes]));
        [y, idx] = ismember(glist, sce.g);
        if ~all(y), error('Runtime error.'); end
        sce.g = sce.g(idx);
        sce.X = sce.X(idx, :);
    end

    % [Xt] = gui.i_transformx(sce.X);
    % if isempty(Xt), return; end

    answer = questdlg('Select method:','','Spearman Correlation','Distance Correlation','Cancel','Spearman Correlation');
    switch answer
        case 'Spearman Correlation'
            fw = gui.gui_waitbar;
            r = corr(t, sce.X.', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
            gui.gui_waitbar(fw);
        case 'Distance Correlation'
            X = log1p(sc_norm(sce.X));
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

    [~, idxp] = maxk(r, 10); 
    selectedg = sce.g(idxp);
    try
        hFig = figure('WindowStyle', 'modal', 'Visible','off');
        hFig.Position(3) = hFig.Position(3) * 1.8;
        gui.i_movegui2parent(hFig, parentfig);

        figure(hFig);
        pkg.i_plot_pseudotimeseries(log2(sce.X+1), ...
            sce.g, t, selectedg);
        hFig.Visible = true;
    catch ME
        if exist('psf1', 'var') && ishandle(hFig)
            close(hFig);
        end
        errordlg(ME.message,'','modal');
    end

end