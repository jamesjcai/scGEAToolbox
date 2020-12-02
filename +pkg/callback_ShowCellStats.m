function callback_ShowCellStats(src,event,sce,h)    
    answer = questdlg('Show cell states?');
    if ~strcmp(answer,'Yes'), return; end    
    
    listitems={'Library Size','Mt-reads Ratio',...
        'Mt-genes Expression','Cell Cycle Phase'};
    for k=1:2:length(sce.list_cell_attributes)
        listitems=[listitems,sce.list_cell_attributes{k}];
    end
    [indx,tf] = listdlg('PromptString',{'Select statistics',...
    '',''},'SelectionMode','single','ListString',listitems);
    if tf==1        
        switch indx
            case 1
                ci=sum(sce.X);
                ttxt="Library Size";
            case 2
                i=startsWith(sce.g,'mt-','IgnoreCase',true);
                lbsz=sum(sce.X,1);
                lbsz_mt=sum(sce.X(i,:),1);
                ci=lbsz_mt./lbsz;
                ttxt="mtDNA%";
            case 3
                idx=startsWith(sce.g,'mt-','IgnoreCase',true);
                n=sum(idx);
                if n>0
                    [ax,bx]=view();
                    if n<=9
                        i_markergenespanel(sce.X,sce.g,sce.s,...
                            sce.g(idx),[],9,ax,bx,'Mt-genes');
                    else
                        i_markergenespanel(sce.X,sce.g,sce.s,...
                            sce.g(idx),[],16,ax,bx,'Mt-genes');
                    end
                else
                    warndlg('No mt-genes found');
                end
                return;
            case 4   % "Cell Cycle Phase";
                if isempty(sce.c_cell_cycle_phase_tx)
                    f = waitbar(0,'Please wait...');
                    pause(.5)
                    waitbar(.67,f,'Processing your data');                
                    [cix]=run_cellcycle(X,genelist);
                    sce.c_cell_cycle_phase_tx=string(cix);
                    waitbar(1,f,'Finishing');
                    pause(1);
                    close(f);
                    
                labels = {'Save cell cycle phase to variable named:'};
                vars = {'c_cell_cycle_phase_tx'};
                values = {sce.c_cell_cycle_phase_tx};
                export2wsdlg(labels,vars,values);
                
                end              
                [ci,tx]=grp2idx(sce.c_cell_cycle_phase_tx);
                ttxt=sprintf('%s|',string(tx));                
            otherwise % other properties
                ttxt=sce.list_cell_attributes{indx-4};
                ci=sce.list_cell_attributes{indx-4+1};                
        end
            [ax,bx]=view();
            delete(h);
            h=i_gscatter3(sce.s,ci,1);
            view(ax,bx);
            title(sce.title);
            if indx==4
                hc=colorbar;
                hc.Label.String=ttxt;
            else
                colorbar off
            end
            colormap default
    end
    
end