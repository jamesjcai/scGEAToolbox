function callback_RunEnrichr_old(src, ~, predefinedlist)

if nargin < 3, predefinedlist = []; end

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    gsorted = natsort(sce.g);

    rng("shuffle");
    n = length(gsorted);
    if isempty(predefinedlist)
        tg = gui.i_inputgenelist(gsorted(randperm(n, min([200, length(gsorted)]))));
    else
        tg = gui.i_inputgenelist(predefinedlist);
    end
    if isempty(tg), return; end

    answer = questdlg('Add background list?','');
    switch answer
        case 'Yes'
            [idx] = gui.i_selmultidlg(sce.g, sce.g, FigureHandle);
            if isempty(idx), return; end
            if idx == 0, return; end
            backgroundlist = sce.g(idx);
        case 'No'
            backgroundlist = [];
        case 'Cancel'
            return;
    end

    if isempty(backgroundlist)
        answer = questdlg('Which Enrichr to use?', '', ...
            'Web-based', ...
            'R/enrichR', ...
            'Python/GSEApy','Web-based');
        
            switch answer
                case 'Web-based'
                    run.web_Enrichr(tg, 200, []);
                    return;
                case 'R/enrichR'
                    [Tbp, Tmf] = run.r_enrichR(tg);
                case 'Python/GSEApy'
                    [Tbp, Tmf] = run.py_GSEApy_enr(tg,[]);
                otherwise
                    return;
            end
    else
        answer = questdlg('Which Enrichr to use?', '', ...
            'Web-based', ...
            'Python/GSEApy','Web-based');
        
            switch answer
                case 'Web-based'
                    run.web_Enrichr_bkg(tg, backgroundlist, 200);
                    return;
                case 'Python/GSEApy'
                    [Tbp, Tmf] = run.py_GSEApy_enr(tg, backgroundlist,...
                        tempdir, false, true);
                otherwise
                    return;
            end
    end
    T = [Tbp; Tmf];
    [~, ~] = gui.i_exporttable(T, true, 'Tenrichrres', 'Enrichr_Results');
    gui.i_viewtable(T, FigureHandle);


    function in_selDataSources
        dsv = ["GO biological process", "GO cellular component", "GO molecular function", ...
               "KEGG", "Reactome", "WikiPathways", "TRANSFAC", "miRTarBase"]; 
        [idx] = gui.i_selmultidlg(dsv, ["GO biological process", "GO molecular function"], FigureHandle);
        
    end

end
