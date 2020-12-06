function callback_GSEA_HVGs(~,~,X,g)
    answer = questdlg('Identify HVGs and then perform GSEA function enrichment analysis?');
    if ~strcmp(answer,'Yes'), return; end 
%     f = waitbar(0,'Please wait...');
%     pause(.5); waitbar(.67,f,'Processing your data');    
t=sc_hvg(X,g);
%     waitbar(1,f,'Finishing');
%     pause(1); close(f);
      

    msgfig=export2wsdlg({'Save HVG table to variable named:'},...
        {'T'},{t});
    uiwait(msgfig)
    
    answer = questdlg('GSEA analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    f = waitbar(0,'Please wait...');
    pause(.5); waitbar(.67,f,'Processing your data');    
tr=run_fgsea(t.genes);
    waitbar(1,f,'Finishing');
    pause(1); close(f);    

        msgfig=export2wsdlg({'Save GSEA table to variable named:'},...
        {'Tr'},{tr});
    uiwait(msgfig)
    
    
        answer = questdlg('GSEA analysis?');
    if ~strcmp(answer,'Yes'), return; end 
    f = waitbar(0,'Please wait...');
    pause(.5); waitbar(.67,f,'Processing your data');    
e_fgseanet(tr);
    waitbar(1,f,'Finishing');
    pause(1); close(f);    
end