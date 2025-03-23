function callback_RunMemento(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    % if ~gui.gui_showrefinfo('Memento [PMID:39454576]', FigureHandle), return; end

    %[wkdir] = gui.i_getwrkdir;
    %if isempty(wkdir), return; end
    extprogname = 'py_memento';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2), return; end
    
    % --------
    a=sprintf('%s vs. %s',cL1{1}, cL2{1});
    b=sprintf('%s vs. %s',cL2{1}, cL1{1});
    answer = gui.myQuestdlg(FigureHandle, 'Which vs. which?','',{a,b},a);
    switch answer
        case a
        case b
            i3=i1; i1=i2; i2=i3;
            cL3=cL1; cL1=cL2; cL2=cL3;
        otherwise
            return;
    end
    % ----------
    X1 = sce.X(:, i1);
    X2 = sce.X(:, i2);
    c = [zeros(size(X1,2),1); ones(size(X2,2),1)];
    scex = SingleCellExperiment([X1 X2], sce.g, [], c);
    scex.c_batch_id = c;
    [succeeded] = run.py_writeh5ad(scex, 'input.h5ad', wkdir);
    if succeeded
        run.py_memento;
    end
