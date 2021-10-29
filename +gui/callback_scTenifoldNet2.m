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
    

answer=questdlg('Use existing networks or reconstruct networks?',...
        'Input Networks','Use existing','Reconstruct','Use existing');
    switch answer
        case 'Use existing'
            a=evalin('base','whos');
            b=struct2cell(a);
            %valididx=ismember(b(4,:),'double');
            %a=a(valididx);
            [indx,tf]=listdlg('PromptString',{'Select two networks:'},...
                'liststring',b(1,:),'SelectionMode','multiple');
            if tf~=1, return; end
            if length(indx)~=2, warndlg('Need two networks.'); return; end
            A0 = evalin('base',a(indx(1)).name);
            A1 = evalin('base',a(indx(2)).name);            
            if ~isequal(size(A0),size(A1))
                errordlg('Two networks should be in the same size.');
                return;
            end            
            [indx,tf]=listdlg('PromptString',{'Select the gene list:'},...
                'liststring',b(1,:),'SelectionMode','single');
            if tf~=1, return; end
            glist = evalin('base',a(indx).name);
            if length(glist)~=size(A0,1)
                errordlg('Not a valid gene list.');
                return;
            end
            fw = gui.gui_waitbar;
            disp('Reading networks...')
            A0sym=0.5*(A0+A0');
            A1sym=0.5*(A1+A1');
            
            disp('Manifold alignment...')
            [aln0,aln1]=i_ma(A0sym,A1sym);
            disp('Differential regulation (DR) detection...')
            T=i_dr(aln0,aln1,glist);
            gui.gui_waitbar(fw);    
            
        case 'Reconstruct'
            [i1,i2]=gui.i_select2grps(sce);
            if length(i1)==1 || length(i2)==1, return; end

            [nsubsmpl,csubsmpl,savegrn]=gui.i_tenifoldnetpara;
            if isempty(nsubsmpl)||isempty(csubsmpl)||isempty(savegrn), return; end
            if csubsmpl>=min([size(sce.X(:,i1),2),size(sce.X(:,i2),2)])
                errordlg('csubsmpl should be smaller than sce.NumCells.');
                return;
            end

            answer123=questdlg('This analysis may take several hours. Continue?');
            if ~strcmpi(answer123,'Yes'), return; end

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

        otherwise
            return;
    end
        

    tstr=matlab.lang.makeValidName(datestr(datetime));
    save(sprintf('T_DRgenes_%s',tstr),'T');
    fprintf('The result has been saved in T_DRgenes_%s.mat\n',tstr);    
    
    figure;
    e_mkqqplot(T);
    gui.i_exporttable(T,true,'T_DRgenes');   
    answer223=questdlg('Run GSEA analysis?');
    if strcmp(answer223,'Yes')
        try
            Tr=ten.e_fgsearun(T);
        catch ME
            errordlg(ME.message);
            return;
        end
        gui.i_exporttable(Tr,true,'T_GSEAres');
        answer323=questdlg('Group GSEA hits?');
        if strcmp(answer323,'Yes')
            ten.e_fgseanet(Tr);
        end
    end
    
    
end
