function callback_GSEA_HVGs(src,~)
    answer = questdlg('Identify HVGs and then perform GSEA function enrichment analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    t=sc_hvg(sce.X,sce.g);
    
    if ~(ismcc || isdeployed)
        msgfig1=export2wsdlg({'Save HVG table to variable named:'},{'T'},{t});
        uiwait(msgfig1)
    else    
        gui.i_exporttable(t,true,'T');
    end

    answer=pkg.timeoutdlg(@(x){questdlg('R/fgeas analysis?')},15);
    % answer = questdlg('GSEA analysis?');
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end 
    
    fw=gui.gui_waitbar;
    tr=run.fgsea(t.genes);
    gui.gui_waitbar(fw);
    
    if ~(ismcc || isdeployed)
        export2wsdlg({'Save GSEA table to variable named:'},{'Tr'},{tr});
        % uiwait(msgfig2)
    else    
        gui.i_exporttable(tr,false,'Tr');
    end
    
    answer=pkg.timeoutdlg(@(x){questdlg('GSEA term network analysis?')},15);
    if strcmp(answer,'No')||strcmp(answer,'Cancel')
        return;
    end
    
    fw=gui.gui_waitbar;
    pkg.e_fgseanet(tr);
    gui.gui_waitbar(fw);
    
    uiwait(helpdlg('Done!',''));
end
