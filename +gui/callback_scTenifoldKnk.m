function callback_scTenifoldKnk(src,~)
    if exist('sctenifoldnet','file')~=2
        errordlg('scTenifoldNet is not installed.');
        disp('To install scTenifoldNet, type:')
        disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
        disp('addpath(''./scTenifoldNet-master/MATLAB'');');
        return;
    end
    if exist('sctenifoldknk','file')~=2
        errordlg('scTenifoldKnk is not installed.');
        disp('To install scTenifoldKnk, type:')
        disp('unzip(''https://github.com/cailab-tamu/scTenifoldKnk/archive/master.zip'');');
        disp('addpath(''./scTenifoldKnk-master/MATLAB'');');
        return;
    end
    
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
            [indx,tf]=listdlg('PromptString',{'Select network varialbe:'},...
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
        answer=questdlg(sprintf('Constructing network and then knocking out gene #%d (%s). Continue?',...
               idx,sce.g(idx)));
    else
        answer=questdlg(sprintf('Knocking out gene #%d (%s) in network (workspace variable %s). Continue?',...
               idx,sce.g(idx),a(indx).name));
    end
    if ~strcmpi(answer,'Yes'), return; end
    
    fw = gui.gui_waitbar;
    try
        if isempty(A0)
            T=sctenifoldknk(sce.X,sce.g,idx,'sorttable',true);
            % T=sortrows(T,'pAdjusted','ascend');
        else            
            A1=A0;
            A1(idx,:)=0;
            [aln0,aln1]=i_ma(A0,A1);
            T=i_dr(aln0,aln1,sce.g,true);
        end
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    gui.i_exporttable(T);
end