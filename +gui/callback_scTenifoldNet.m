function callback_scTenifoldNet(src,~)

    if exist('sctenifoldnet','file')~=2
        errordlg('scTenifoldNet is not installed.');
        return;
    end
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [i1,i2]=gui.i_select2grps(sce);
    if length(i1)==1 || length(i2)==1, return; end

    answer=questdlg('This analysis may take several hours. Continue?');
    if ~strcmpi(answer,'Yes'), return; end
    savegrn=false;
    answer=questdlg('Save constructed networks (to current folder)?');
    switch answer
        case 'Yes'
            savegrn=true;
        case 'No'
            savegrn=false;
        otherwise
            return;
    end    
    fw = gui.gui_waitbar;
    try
    [T]=sctenifoldnet(sce.X(:,i1),sce.X(:,i2),sce.g,...
                    'nsubsmpl',10,'savegrn',savegrn);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    figure;
    e_mkqqplot(T);    
    gui.i_exporttable(T,true);
end