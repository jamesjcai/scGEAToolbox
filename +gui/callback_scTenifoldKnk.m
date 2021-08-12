function callback_scTenifoldKnk(src,~)

%     if exist('sctenifoldnet','file')~=2
%         errordlg('scTenifoldNet is not installed.');
%         disp('To install scTenifoldNet, type:')
%         disp('unzip(''https://github.com/cailab-tamu/scTenifoldNet/archive/master.zip'');');
%         disp('addpath(''./scTenifoldNet-master/MATLAB'');');
%         return;
%     end
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    
    a=evalin('base','whos');
    b=struct2cell(a);
    %valididx=ismember(b(4,:),'double');
    %a=a(valididx);
    [indx,tf]=listdlg('PromptString',{'Select network varialbe:'},...
        'liststring',b(1,:),'SelectionMode','single');
    if tf==1
        A = evalin('base',a(indx).name);
    else
        return;
    end    

    [m,n]=size(A);
    if m~=n || n~=length(sce.g), errordlg('Not a valid network.'); return; end    
    
    gsorted=sort(sce.g);
    [indx2,tf] = listdlg('PromptString',{'Select the KO gene:'},...
        'SelectionMode','single','ListString',gsorted);
    if tf==1
        [~,idx]=ismember(gsorted(indx2),sce.g);
    else
        return;
    end

    answer=questdlg(sprintf('Knocking out gene #%d (%s) in network (workspace variable %s). Continue?',...
           idx,sce.g(idx),a(indx).name));
    if ~strcmpi(answer,'Yes'), return; end
    
    fw = gui.gui_waitbar;
    try
    % [T]=sctenifoldknk(A,idx);
    catch ME
        gui.gui_waitbar(fw);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    % gui.i_exporttable(T);
end