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

    [nsubsmpl,csubsmpl,savegrn]=gui.i_tenifoldnetpara;
    if isempty(nsubsmpl)||isempty(csubsmpl)||isempty(savegrn), return; end
    if csubsmpl>=min([size(sce.X(:,i1),2),size(sce.X(:,i2),2)])
        errordlg('csubsmpl should be smaller than sce.NumCells.');
        return;
    end
    
    answer=questdlg('This analysis may take several hours. Continue?');
    if ~strcmpi(answer,'Yes'), return; end
    
    fw = gui.gui_waitbar;
    try
        fprintf('\n');
        disp('[T]=ten.sctenifoldnet(X1,X2,g,''nsubsmpl'',10,''csubsmpl'',500,''savegrn'',true);')
        disp('[T]=ten.sctenifoldnet(sce.X(:,idx1),sce.X(:,idx2),sce.g,''nsubsmpl'',10,''csubsmpl'',500,''savegrn'',true);')
        [T]=ten.sctenifoldnet(sce.X(:,i1),sce.X(:,i2),sce.g,...
           'nsubsmpl',nsubsmpl,'csubsmpl',csubsmpl,'savegrn',savegrn);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    tstr=matlab.lang.makeValidName(datestr(datetime));
    save(sprintf('T_DRgenes_%s',tstr),'T');
    fprintf('The result has been saved in T_DRgenes_%s.mat\n',tstr);    
    gui.gui_waitbar(fw);
    figure;
    e_mkqqplot(T);
    gui.i_exporttable(T,true,'T_DRgenes');   
    answer=questdlg('Run GSEA analysis?');
    if strcmp(answer,'Yes')
        try
            Tr=ten.e_fgsearun(T);
        catch ME
            errordlg(ME.message);
            return;
        end
        gui.i_exporttable(Tr,true,'T_GSEAres');
        answer2=questdlg('Group GSEA hits?');
        if strcmp(answer2,'Yes')
            ten.e_fgseanet(Tr);
        end
    end
end
