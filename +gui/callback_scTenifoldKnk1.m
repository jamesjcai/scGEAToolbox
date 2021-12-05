function callback_scTenifoldKnk1(src,~)
    import ten.*
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    answer=questdlg('Use existing network or reconstruct network?',...
        'Input Network','Use existing','Reconstruct','Use existing');
    switch answer
        case 'Use existing'
            a=evalin('base','whos');
            b=struct2cell(a);
            valididx=false(length(a),1);
            for k=1:length(a)
                if max(a(k).size)==sce.NumGenes && min(a(k).size)==sce.NumGenes
                    valididx(k)=true;
                end
            end
            if ~any(valididx)
                warndlg('Workspace contains no network varible.');
                return;
            else
                %valididx=ismember(b(4,:),'double');
                a=a(valididx);
                b=b(:,valididx);
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
            end
        case 'Reconstruct'
            try
                ten.check_tensor_toolbox;
            catch ME        
                errordlg(ME.message);
                return;
            end            
            A0=[];
            uiwait(helpdlg("Network will be constructed. Now select the gene to be knocked out..."));
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
    

    
        if isempty(A0)
            try
                fw = gui.gui_waitbar;
                [T,A0]=ten.sctenifoldknk(sce.X,sce.g,idx,'sorttable',true);                
                gui.gui_waitbar(fw);
             catch ME
                gui.gui_waitbar(fw);
                errordlg(ME.message);
                return;
            end
            isreconstructed=true;
        else
            doit=false;
            if sum(A0(idx,:)~=0)==0
                s=sprintf('KO gene (%s) has no link or too few links (n<50) with other genes.',...
                          sce.g(idx));
                warndlg(s);
                return;
            elseif sum(A0(idx,:)~=0)<50
                s=sprintf('KO gene (%s) has too few links (n=%d) with other genes. Continue?',...
                          sce.g(idx),sum(A0(idx,:)~=0));
                answer11 = questdlg(s);
                switch answer11
                    case 'Yes'
                        doit=true;                        
                    case 'No'
                        return;
                    case 'Cancel'
                        return;
                    otherwise
                        return;
                end
            else
                doit=true;                
            end
            
            if doit
                try
                    fw = gui.gui_waitbar;
                    disp('>> [T]=ten.i_knk(A0,targetgene,genelist,true);')
                    [T]=ten.i_knk(A0,idx,sce.g,true);
                    gui.gui_waitbar(fw);
                catch ME
                    gui.gui_waitbar(fw);
                    errordlg(ME.message);
                    return;
                end
            end
            %A1=A0;
            %A1(idx,:)=0;
            %[aln0,aln1]=i_ma(A0,A1);
            %T=i_dr(aln0,aln1,sce.g,true);
            isreconstructed=false;
        end    
    
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
    disp('Tn=ten.e_fgseanet(Tf);');
    disp('===============================');
end
