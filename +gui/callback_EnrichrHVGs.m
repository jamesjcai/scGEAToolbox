function callback_EnrichrHVGs(src,~)
    answer = questdlg('Identify HVGs and perform function enrichment analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
  
    answer = questdlg('Which method?','Select Method', ...        
        'Brennecke et al. (2013)','Splinefit Method',...
        'Brennecke et al. (2013)');
    switch answer
        case 'Brennecke et al. (2013)'
            fw = gui.gui_waitbar;
            t=sc_hvg(sce.X,sce.g,true,true);
            gui.gui_waitbar(fw);

            if ~(ismcc || isdeployed)
                msgfig1=export2wsdlg({'Save HVG table to variable named:'},{'T'},{t});
                uiwait(msgfig1)
            else
                gui.i_exporttable(t,true,'T');
            end
        
            gui.i_enrichtest(t.genes);

        case 'Splinefit Method'
            fw = gui.gui_waitbar;
            try
                gui.sc_scatter3genes(sce.X,sce.g);
            catch ME
                gui.gui_waitbar(fw,true);
                errordlg(ME.message);
            end
            gui.gui_waitbar(fw,true);
        otherwise
            return;
    end
end
