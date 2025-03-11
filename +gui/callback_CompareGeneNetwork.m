function callback_CompareGeneNetwork(src, ~)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    
    [i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
    if isscalar(i1) || isscalar(i2)
        if i1 == 0 || i2 == 0, return; end
    end
    
    [glist] = gui.i_selectngenes(sce,[],FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(glist, sce.g);
    if ~all(y), error('Selected gene(s) not in the gene list of data.'); end
    fprintf("%s\n", glist)
    
    switch gui.myQuestdlg(FigureHandle, "Select algorithm:",'',...
            {'PC Regression','Chaterjee Correlation'},'PC Regression')
        case 'PC Regression'
            [Xt] = gui.i_transformx(sce.X, true, 5, FigureHandle);
            if isempty(Xt), return; end
            
            x1 = Xt(i, i1);
            x2 = Xt(i, i2);
            
            fw = gui.gui_waitbar;
            A1 = sc_pcnet(x1, 3, false, true, false);
            A2 = sc_pcnet(x2, 3, false, true, false);
            gui.gui_waitbar(fw);
        case 'Chaterjee Correlation'
            [Xt] = gui.i_transformx(sce.X, true, 3, FigureHandle);
            if isempty(Xt), return; end
            
            x1 = Xt(i, i1);
            x2 = Xt(i, i2);

            fw = gui.gui_waitbar_adv;
            n1 = size(x1,1);
            A1 = zeros(n1);
            n2 = size(x2,1);
            A2 = zeros(n2);
            for k=1:n1
                gui.gui_waitbar_adv(fw,k/n1);
                for l=1:n1
                    if k~=l
                        A1(k,l) = pkg.e_xicor(x1(k,:),x1(l,:));
                    end
                end
                for l=1:n2
                    if k~=l
                        A2(k,l) = pkg.e_xicor(x2(k,:),x2(l,:));
                    end
                end
            end
            gui.gui_waitbar_adv(fw);
        otherwise
            return;
    end
    pause(1)
    stitle = sprintf('%s vs. %s', cL1{1}, cL2{1});
    try
        sc_grnview2(A1, A2, glist, stitle, FigureHandle);
    catch ME
        errordlg(ME.message);
    end
end
