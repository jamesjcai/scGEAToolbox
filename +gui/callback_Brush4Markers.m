function callback_Brush4Markers(src, event)


    [FigureHandle, sce] = gui.gui_getfigsce(src);

    if ~gui.i_installed('stats'), return; end

    switch gui.myQuestdlg(FigureHandle, 'Select method:','',...
            {'Lasso Regression','Logistic Regression 🐢 '}, ...
             'Lasso Regression')
        case 'Lasso Regression'
            uselasso=true;
        case 'Logistic Regression 🐢 '
            uselasso=false;
        otherwise
            return;
    end
    i_Brush4MarkersLASSO(src, event, sce, uselasso);
end


function i_Brush4MarkersLASSO(src, ~, sce, uselasso)
    [FigureHandle] = gui.gui_getfigsce(src);

    if nargin < 3
        uselasso = true; 
    end
    
    axesh = FigureHandle.findobj('type', 'Axes');
    [axx, bxx] = view(axesh);
    
    % [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    
    assert(isequal(axesh.findobj('type', 'Scatter'), ...
        FigureHandle.findobj('type', 'Scatter')))
    %axesh.Children(1)
    %isequal(axesh.findobj('type','Scatter'),axesh.Children(2))
    h = axesh.findobj('type', 'Scatter');
    ptsSelected = logical(h.BrushData.');
    
    
    if ~any(ptsSelected)
        answer=gui.myQuestdlg(FigureHandle, 'No cells are brushed/selected. You can select cells by a grouping variable. Continue?','');
        if ~strcmp(answer,'Yes'), return; end
        [ptsSelected] = gui.i_select1classcells(sce, false);
        if isempty(ptsSelected), return; end
        if all(ptsSelected)
            gui.myWarndlg(FigureHandle, "All cells are in the same group.");
            return;
        end
    else
        % assignin('base', 'ptsSelected', ptsSelected);
        [ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce);
        if ~letdoit, return; end
    end
    
    
    [numfig] = gui.i_inputnumg(500);
    if isempty(numfig), return; end
    
    
    if uselasso, fw = gui.myWaitbar(FigureHandle); end
    y = double(ptsSelected);
    sce.c = 1 + ptsSelected;
    X = sce.X';
    
    % uselasso = true;
    
    try
        if issparse(X), X = full(X); end
    %    assignin("base","X",X);
    %    assignin("base","y",y);
        if uselasso
            [B] = lasso(X, y, 'DFmax', numfig*3, 'MaxIter', 1e3);
            [~, ix] = min(abs(sum(B > 0)-numfig));
            b = B(:, ix);
            idx = b > 0;
        else
            %mdl = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');
            %B = mdl.Coefficients.Estimate;
            %[~, idx] = mink(B, numfig);
            idx = LRDETest(X, y, numfig);
        end
    catch ME
        if uselasso, gui.myWaitbar(FigureHandle, fw, true); end
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        rethrow(ME);
    end
    
    if ~any(idx)
       if uselasso, gui.myWaitbar(FigureHandle, fw); end
        gui.myWarndlg(FigureHandle, 'No marker found');
        return;
    end
    
    if uselasso
        markerlist = sce.g(idx);
        [~, jx] = sort(b(idx), 'descend');
        markerlist = markerlist(jx);
    else
        markerlist = sce.g(idx);
    end
    
    fprintf('%d marker genes: ', length(markerlist));
    fprintf('%s ', markerlist)
    fprintf('\n')

    n = length(markerlist);
    y = cell(n,1);
    for k=1:n
        y{k} = sce.X(sce.g == markerlist(k), :);
    end

    gui.sc_uitabgrpfig_expplot(y, markerlist, sce.s, FigureHandle, [axx, bxx]);
    if uselasso, gui.myWaitbar(FigureHandle, fw); end
end

function idx = LRDETest(X, y, k)
    n = size(X, 1);
    p_val = zeros(n, 1);
    fw = gui.myWaitbar(FigureHandle);
    % Calculate p-values for each row of data_use
    for x = 1:size(X, 1)
        if mod(x,5)==0
            gui.myWaitbar(FigureHandle, fw, false, '', '', x/n);
        end
        model_data = table(X(:,x), y(:), 'VariableNames', {'GENE', 'Group'});
        fmla = 'Group ~ GENE';
        fmla2 = 'Group ~ 1';
        % Fit models and compute likelihood ratio test
        model1 = fitglm(model_data, fmla, 'Distribution', 'binomial');
        model2 = fitglm(model_data, fmla2, 'Distribution', 'binomial');
        lrtest_stat = 2 * (model1.LogLikelihood - model2.LogLikelihood);
        p_val(x) = chi2cdf(lrtest_stat, model1.NumPredictors - model2.NumPredictors, 'upper');
    end
    [~, idx] = mink(p_val, k);
    gui.myWaitbar(FigureHandle, fw);
end
