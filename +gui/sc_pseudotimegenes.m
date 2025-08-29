function sc_pseudotimegenes(sce, t, parentfig)

if nargin<3, parentfig=[]; end
t = t(:);

    [K, usehvgs] = gui.i_gethvgnum(sce, parentfig);
    
    if usehvgs
        T = sc_splinefit(sce.X, sce.g);
        glist = T.genes(1:min([K, sce.NumGenes]));
        [y, idx] = ismember(glist, sce.g);
        if ~all(y), error('Runtime error.'); end        
        sce.X = sce.X(idx, :);
        sce.g = sce.g(idx);
    end

    % [Xt] = gui.i_transformx(sce.X, [], [], parentfig);
    % if isempty(Xt), return; end
    X = sce.X;
    try
        if issparse(X), X = full(X); end
    catch ME
        disp(ME.message);
        disp('Keep using sparse X.');
    end

    

    answer = gui.myQuestdlg(parentfig, 'Select method:','',{'Spearman Correlation', ...
        'Distance Correlation','Cancel'},'Spearman Correlation');
    switch answer
        case 'Spearman Correlation'
            fw = gui.myWaitbar(parentfig);
            r = corr(t, X.', 'type', 'spearman'); % Calculate linear correlation between gene expression profile and T
            gui.myWaitbar(parentfig, fw);
        case 'Distance Correlation'
            Xn = log1p(sc_norm(X));
            r = zeros(size(X,1),1);
            n = length(r);
            if n > 100
                fw = gui.myWaitbar(parentfig);
                gui.myWaitbar(parentfig, fw, false, '', '', fw,100/n);
                for k=1:length(r)
                    if mod(k,100) == 0
                        gui.myWaitbar(parentfig, fw, false, '', '', k/length(r)); 
                    end
                    r(k) = pkg.distcorr(Xn(k,:).', t(:)); % Calculate linear correlation between gene expression profile and T
                end
                gui.myWaitbar(parentfig, fw);
            else
                fw = gui.myWaitbar(parentfig);
                for k=1:length(r)
                    r(k) = pkg.distcorr(Xn(k,:).', t(:)); % Calculate linear correlation between gene expression profile and T
                end
                gui.myWaitbar(parentfig, fw);
            end
        case 'Cancel'
            return; 
        otherwise
            return;
    end

    [~, idxp] = maxk(r, 10); 
    selectedg = sce.g(idxp);
    try
        hx=gui.myFigure(parentfig);
        hFig=hx.FigHandle;
        hFig.Position(3) = hFig.Position(3) * 1.8;
        hx.show(parentfig);
        figure(hFig);
        pkg.i_plot_pseudotimeseries(log1p(X), ...
            sce.g, t, selectedg);        
    catch ME
        if exist('psf1', 'var') && ishandle(hFig)
            close(hFig);
        end
        gui.myErrordlg(parentfig, ME.message,'','modal');
    end

end