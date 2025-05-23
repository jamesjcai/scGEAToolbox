function callback_scTenifoldNet2lite(src, ~)

import ten.*
import pkg.*

[FigureHandle, sce] = gui.gui_getfigsce(src);
    
if ~gui.gui_showrefinfo('scTenifoldNet [PMID:33336197]', FigureHandle), return; end

extprogname = 'ml_scTenifoldNet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end


[i1, i2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end

    fw = gui.myWaitbar(FigureHandle);
    disp('Constructing networks (1/2) ...')
    X = sc_norm(sce.X);
    X = log1p(X);
    X0 = X(:, i1);
    X1 = X(:, i2);
    A0 = sc_pcnetpar(X0);
    disp('Constructing networks (2/2) ...')
    A1 = sc_pcnetpar(X1);
    A0sym = 0.5 * (A0 + A0');
    A1sym = 0.5 * (A1 + A1');
    
    disp('Manifold alignment...')
    [aln0, aln1] = i_ma(A0sym, A1sym);
    disp('Differential regulation (DR) detection...')
    glist = sce.g;
    T = i_dr(aln0, aln1, glist);
    gui.myWaitbar(FigureHandle, fw);
    
    tstr = matlab.lang.makeValidName(string(datetime));
    b = 'sctenifoldnet_outs';
    a = sprintf('output_%s', tstr);
    if ~exist(fullfile(wkdir, b), 'dir')
        mkdir(fullfile(wkdir, b));
        pause(1);
    end
    f1 = fullfile(wkdir, b, a);
    save(f1, 'T', 'A0', 'A1', 'glist');
    a = sprintf('output_%s.xlsx', tstr);
    f1 = fullfile(wkdir, b, a);
    writetable(T, f1, 'FileType', 'spreadsheet');
    fprintf('The result has been saved in %s\n', f1);
    bx = sprintf('The result has been saved in %s. Open the folder to locate it?', f1);
    if strcmp('Yes', gui.myQuestdlg(FigureHandle, bx))
            winopen(fullfile(wkdir, b));
    end

end
