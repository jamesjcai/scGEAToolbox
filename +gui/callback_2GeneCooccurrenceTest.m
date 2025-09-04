function callback_2GeneCooccurrenceTest(src, ~)

    [FigureHandle, sce_ori] = gui.gui_getfigsce(src);
    sce = copy(sce_ori);


    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    if isscalar(glist)
        gui.myErrordlg(FigureHandle,'Two genes are needed for the anlaysis.');
        return;
    end
    if length(glist) == 2
        g1 = glist(1);
        g2 = glist(2);
    elseif length(glist)>2
        g1 = glist(1);
        g2 = glist(2);        
        if ~strcmp('Yes', gui.myQuestdlg(FigureHandle, sprintf(['Only two ' ...
            'selected genes, %s and %s, will be used for the analysis. Continue?'], g1, g2),'',[],[],'warning'))
            return;
        end
    end



    answer = gui.myQuestdlg(FigureHandle, 'All cells, or select cells?', ...
        'Apply Analysis to All Cells?',{'All Cells','Select Cells'});
    switch answer
        case 'All Cells'
            fw = gui.myWaitbar(FigureHandle);
        case 'Select Cells'
            [thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
            if isempty(thisc), return; end
            [~, cL] = findgroups(string(thisc));
            [idx] = gui.i_selmultidialog(cL, [], FigureHandle);
            if isempty(idx), return; end
            fw = gui.myWaitbar(FigureHandle);
            sce = sce.selectcells(ismember(thisc, cL(idx)));  %#OK
        otherwise
            return;
    end


    a=0+sce.X(sce.g==g1,:)>0;
    b=0+sce.X(sce.g==g2,:)>0;
    
    % [chi2_stat, p_value, contingency_table, expected_freq] = pkg.ai_chi2binarytest(a,b);
    % Capture the printed output into a variable
    
    try
        pkg.ai_chi2binarytest
    catch
    end
   
    capturedText = evalc('pkg.ai_chi2binarytest(a,b);');
    
    Xm = sc_impute(sce.X);


    gui.myWaitbar(FigureHandle, fw);

    figure;
    scatter(Xm(sce.g==g1,:),Xm(sce.g==g2,:));
    xlabel(g1);
    ylabel(g2);    

    %if gui.i_isuifig(FigureHandle)
        %gui.myInputdlg({'Test Info:'}, 'Output Viewer', {char(capturedText)}, FigureHandle);
    %    gui.myTextareadlg(FigureHandle, {'Test Info:'}, 'Output Viewer', {capturedText}, true);
    %else
        inputdlg('Test Info:', 'Output Viewer', [15, 80], {char(capturedText)});
    %end   

end
