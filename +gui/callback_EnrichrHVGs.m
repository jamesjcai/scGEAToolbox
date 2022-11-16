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

%             answer1=pkg.timeoutdlg(@(x){questdlg('Which functional enrichment analysis do you want to use?','Analysis Method', ...
%                 'Enrichr','GOrilla','Enrichr+GOrilla','Enrichr')},15);
%             if isempty(answer1), return; end
%             switch answer1
%                 case 'Enrichr'
%                     run.Enrichr(t.genes,500);
%                 case 'GOrilla'
%                     run.GOrilla(t.genes(1:500));
%                 case 'Enrichr+GOrilla'
%                     run.Enrichr(t.genes,500);
%                     run.GOrilla(t.genes(1:500));
%                 otherwise
%                     return;
%             end
            
            
%             fw=gui.gui_waitbar;
%             tr=run.fgsea(t.genes);
%             gui.gui_waitbar(fw);
%             
%             if ~(ismcc || isdeployed)
%                 export2wsdlg({'Save GSEA table to variable named:'},{'Tr'},{tr});
%                 % uiwait(msgfig2)
%             else    
%                 gui.i_exporttable(tr,false,'Tr');
%             end            
%            answer=pkg.timeoutdlg(@(x){questdlg('GSEA term network analysis?')},15);
%            if strcmp(answer,'No')||strcmp(answer,'Cancel')
%                return;
%            end
%            fw=gui.gui_waitbar;
%            pkg.e_fgseanet(tr);
%            gui.gui_waitbar(fw);
%            uiwait(helpdlg('Done!',''));

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
