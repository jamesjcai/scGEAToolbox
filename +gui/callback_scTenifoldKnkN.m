function callback_scTenifoldKnkN(src,~)
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
    
%     gsorted=sort(sce.g);
%     [indx2,tf] = listdlg('PromptString',{'Select the KO gene:'},...
%         'SelectionMode','single','ListString',gsorted);
%     if tf==1
%         [~,idx]=ismember(gsorted(indx2),sce.g);
%     else
%         return;
%     end
    
%     if isempty(A0)
%         answer=questdlg(sprintf('Ready to construct network and then knock out gene #%d (%s). Continue?',...
%                idx,sce.g(idx)));
%     else
%         answer=questdlg(sprintf('Ready to knock out gene #%d (%s) from network (%s). Continue?',...
%                idx,sce.g(idx),a(indx).name));
%     end
%     
%     if ~strcmpi(answer,'Yes'), return; end

        if isempty(A0)
            try
                fw = gui.gui_waitbar;
                [A0]=sc_pcnetdenoised(sce.X);
                [F]=ten.knk3_buildPerturbationLandscape(A0,sce.g);
                gui.gui_waitbar(fw);
             catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return;
            end
            isreconstructed=true;
        else
%                try
                    fw = gui.gui_waitbar;                    
                    [F]=ten.knk3_buildPerturbationLandscape(A0,sce.g);
                    gui.gui_waitbar(fw);
%                 catch ME
%                     gui.gui_waitbar(fw);
%                     errordlg(ME.message);
%                     return;
%                 end
            %A1=A0;
            %A1(idx,:)=0;
            %[aln0,aln1]=i_ma(A0,A1);
            %T=i_dr(aln0,aln1,sce.g,true);
            isreconstructed=false;
        end
        
    if isreconstructed
        labels = {'Save network to variable named:',...
            'Save perturbation score to variable named:',...
            'Save gene list to variable named:'}; 

        vars = {'A0','F','g'};
        values = {A0,F,sce.g};
    else
        labels = {'Save perturbation score to variable named:',...
            'Save gene list to variable named:'}; 
        vars = {'F','g'};
        values = {F,sce.g};
    end
    waitfor(export2wsdlg(labels,vars,values));
end
