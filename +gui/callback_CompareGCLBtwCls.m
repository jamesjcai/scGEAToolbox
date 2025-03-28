function callback_CompareGCLBtwCls(src, ~)
% callback_CompareGCLBtwCls compares Global Coordination Level (GCL) in genes to show differences between cell groups.
%
% Inputs:
%   src - Source input for the GUI.

%gui.gui_showrefinfo('GCL Analysis [PMID:33139959]', FigureHandle);

    [FigureHandle, sce] = gui.gui_getfigsce(src);
        
    if ~gui.gui_showrefinfo('GCL Analysis [PMID:33139959]', FigureHandle), return; end
    answer = gui.myQuestdlg(FigureHandle, 'This function compares GCL of genes to show differences between cell groups. Continue?', '');
    if ~strcmp(answer, 'Yes'), return; end
      
    [thisc, clabel] = gui.i_select1class(sce,[],[],[],FigureHandle);
    if isempty(thisc), return; end
            
    n = length(unique(thisc));
    if n == 1
        answer = gui.myQuestdlg(FigureHandle, "All cells are in the same group. Continue?", "");
        switch answer
            case 'Yes'
            otherwise
                return;
        end
    end
    
    %'Global Coordination Level (GCL) [PMID:33139959]'
    [c, cL] = grp2idx(thisc);
    
    N = 50;
    
    t = zeros(N*n, 1);
    V = zeros(N*n, 1);
    fw = gui.myWaitbar(FigureHandle);
    
    for k = 1:n
        fprintf('Working on %s: %s ... %d of %d\n', clabel, cL{k}, k, n);
        idx = (k - 1) * N + 1:k * N;
        [v] = run.ml_GCL(sce.X(:, k == c), N);
        t(idx, :) = k;
        V(idx, :) = v;
    end
    gui.myWaitbar(FigureHandle, fw);
    
    y = V;
    thisc = cL(t);
    ttxt = 'GCL';
    gui.i_violinplot(y, thisc, ttxt, true, [], [], FigureHandle);


end