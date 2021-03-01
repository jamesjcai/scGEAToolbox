function callback_GSEA_HVGs(src,~)
    answer = questdlg('Identify HVGs and then perform GSEA function enrichment analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    t=sc_hvg(sce.X,sce.g);
    msgfig1=export2wsdlg({'Save HVG table to variable named:'},...
        {'T'},{t});
    % uiwait(msgfig)
    
    answer=pkg.timeoutdlg(@(x){questdlg('GSEA analysis?')},15);
    % answer = questdlg('GSEA analysis?');
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end 
    
    fw=gui.gui_waitbar;
    tr=run_fgsea(t.genes);
    gui.gui_waitbar(fw);
    msgfig2=export2wsdlg({'Save GSEA table to variable named:'},...
        {'Tr'},{tr});
    % uiwait(msgfig)
    
    
    answer=pkg.timeoutdlg(@(x){questdlg('GSEA term network analysis?')},15);
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end
    
    fw=gui.gui_waitbar;
    e_fgseanet(tr);
    gui.gui_waitbar(fw);
    
            a=helpdlg({'Done!'});
             uiwait(a)
end
