function [cs] = e_cellscore(sce, posg, methodid, showwaitbar, FigureHandle)

if nargin < 5, FigureHandle = []; end

if nargin < 4, showwaitbar=true; end

cs = [];
if nargin < 3 || isempty(methodid)
    [~, methodid] = gui.i_pickscoremethod(0, FigureHandle);
    % answer = gui.myQuestdlg(FigureHandle, 'Select algorithm:', ...
    %     'Select Method', ...
    %     'UCell [PMID:34285779]', 'AddModuleScore/Seurat', ...
    %     'UCell [PMID:34285779]');
    % switch answer
    %     case 'AddModuleScore/Seurat'
    %         methodid = 2;
    %     case 'UCell [PMID:34285779]'
    %         methodid = 1;
    %     otherwise
    %         return;
    % end
    if isempty(methodid), return; end
end

if showwaitbar, fw = gui.gui_waitbar; end
try
    if methodid == 1
        [cs] = sc_cellscore_ucell(sce.X, sce.g, posg);
    elseif methodid == 2
        [cs] = sc_cellscore_admdl(sce.X, sce.g, posg);
    end
catch ME
    if showwaitbar, gui.gui_waitbar(fw, true); end
    errordlg(ME.message);
    return;
end
if showwaitbar, gui.gui_waitbar(fw); end

if showwaitbar
    posg = sort(posg);
    isexpressed = ismember(upper(posg), upper(sce.g));
    fprintf('\n=============\n%s\n-------------\n', 'Genes');
    for k = 1:length(posg)
        if isexpressed(k)
            fprintf('%s*, ', posg(k));
        else
            fprintf('%s, ', posg(k));
        end
        if mod(k, 10) == 0 || k == length(posg)
            fprintf('\n');
        end
    end
    fprintf('=============\n*Expressed genes (n = %d)\n', ...
        sum(isexpressed));
end

end