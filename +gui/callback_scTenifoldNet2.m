function callback_scTenifoldNet2(src,~)
    import ten.*
    import pkg.*
    try
        ten.check_tensor_toolbox;
    catch ME
        errordlg(ME.message);
        return;
    end
    
%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [i1,i2]=gui.i_select2grps(sce);
    if length(i1)==1 || length(i2)==1, return; end

    answer=questdlg('This analysis may take several hours. Continue?');
    if ~strcmpi(answer,'Yes'), return; end

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
    [T]=ten.sctenifoldnet(sce.X(:,i1),sce.X(:,i2),sce.g,...
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
    answer=questdlg('Run GSEA analysis?');
    if strcmp(answer,'Yes')
        try
            Tr=ten.e_fgsearun(T);
        catch ME
            errordlg(ME.message);
            return;
        end
        gui.i_exporttable(Tr,true);
        answer2=questdlg('Group GSEA hits?');
        if strcmp(answer2,'Yes')
            ten.e_fgseanet(Tr);
        end
    end
end
