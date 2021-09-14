function callback_scTenifoldKnk(src,~)
    import ten.*
%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end
%     if exist('sctenifoldknk','file')~=2
%         errordlg('scTenifoldKnk is not installed.');
%         disp('To install scTenifoldKnk, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldKnk/archive/master.zip'');');
%         disp('addpath(''./scTenifoldKnk-master/MATLAB'');');
%         return;
%     end
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    answer=questdlg('Use existing network or reconstruct network?',...
        'Input Network','Use existing','Reconstruct','Use existing');
    switch answer
        case 'Use existing'
            a=evalin('base','whos');
            b=struct2cell(a);
            %valididx=ismember(b(4,:),'double');
            %a=a(valididx);
            [indx,tf]=listdlg('PromptString',{'Select network variable:'},...
                'liststring',b(1,:),'SelectionMode','single');
            if tf==1
                A0 = evalin('base',a(indx).name);
            else
                return;
            end
            [m,n]=size(A0);
            if m~=n || n~=length(sce.g)
                errordlg('Not a valid network.'); 
                return;
            end    
        case 'Reconstruct'
            try
                ten.check_tensor_toolbox;
            catch ME        
                errordlg(ME.message);
                return;
            end            
            A0=[];
        otherwise
            return;
    end
    
    gsorted=sort(sce.g);
    [indx2,tf] = listdlg('PromptString',{'Select the KO gene:'},...
        'SelectionMode','single','ListString',gsorted);
    if tf==1
        [~,idx]=ismember(gsorted(indx2),sce.g);
    else
        return;
    end
    
    if isempty(A0)
        answer=questdlg(sprintf('Ready to construct network and then knock out gene #%d (%s). Continue?',...
               idx,sce.g(idx)));
    else
        answer=questdlg(sprintf('Ready to knock out gene #%d (%s) from network (%s). Continue?',...
               idx,sce.g(idx),a(indx).name));
    end
    
    if ~strcmpi(answer,'Yes'), return; end
    
    fw = gui.gui_waitbar;
    try
        if isempty(A0)
            [T,A0]=ten.sctenifoldknk(sce.X,sce.g,idx,'sorttable',true);
            % T=sortrows(T,'pAdjusted','ascend');
            isreconstructed=true;
        else
            [T]=i_knk(A0,idx,sce.g);
            %A1=A0;
            %A1(idx,:)=0;
            %[aln0,aln1]=i_ma(A0,A1);
            %T=i_dr(aln0,aln1,sce.g,true);
            isreconstructed=false;
        end
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    
    if isreconstructed
        labels = {'Save network to variable named:'}; 
        vars = {'A0'};
        values = {A0};        
        waitfor(export2wsdlg(labels,vars,values));
    end
    
    gui.i_exporttable(T);
    disp('Downstream Analysis Options:');
    disp('===============================');
    disp('run.Enrichr(T.genelist(1:200));');
    disp('Tf=ten.e_fgsearun(T);');
    disp('Tn=e_fgseanet(Tf);');
    disp('===============================');
end
