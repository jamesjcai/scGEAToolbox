function callback_GSEA_HVGs(~,~,X,g)
    answer = questdlg('Identify HVGs and then perform GSEA function enrichment analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    t=sc_hvg(X,g);
    msgfig=export2wsdlg({'Save HVG table to variable named:'},...
        {'T'},{t});
    uiwait(msgfig)
    
    answer=pkg.timeoutdlg(@(x){questdlg('GSEA analysis?')},15);
    % answer = questdlg('GSEA analysis?');
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end 
    
    fw=pkg.gui_waitbar;
    tr=run_fgsea(t.genes);
    pkg.gui_waitbar(fw);
    msgfig=export2wsdlg({'Save GSEA table to variable named:'},...
        {'Tr'},{tr});
    uiwait(msgfig)
    
    
    answer=pkg.timeoutdlg(@(x){questdlg('GSEA term network analysis?')},15);
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end 
    
    fw=pkg.gui_waitbar;
    e_fgseanet(tr);
    pkg.gui_waitbar(fw);
end