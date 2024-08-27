function [isok] = i_enrichrprep
    isok = false;
    %{
    runenrichr = questdlg('Run Enrichr (R required) with top 250 DE genes? Results will be saved in the output Excel files.','');
    if strcmp(runenrichr,'Cancel'), return; end
    isok = false;
    if strcmp(runenrichr, 'Yes')
        [isok] = gui.i_commoncheck_R('r_enrichR');
        if ~isok
            answer = questdlg('R Environment Error: It seems that your R environment is not configured correctly. This means that Gene Function Enrichment Analysis cannot be performed. Continue?',''); 
            if ~strcmp(answer,'Yes'), return; end
        end
    end
    %}
    
    runenrichr = questdlg('Run Enrichr (Python required) with top 250 DE genes? Results will be saved in the output Excel files.','');
    if strcmp(runenrichr,'Cancel'), return; end
    if strcmp(runenrichr, 'Yes')
        [good] = gui.i_commoncheck_Python('py_GSEApy_enr');
        if good
            isok = true;
        else
            answer = questdlg('Python Environment Error: It seems that your Python environment is not configured correctly. This means that Gene Function Enrichment Analysis cannot be performed. Continue?',''); 
            if strcmp(answer,'Yes')
                isok = [];    % continue without enrichr
            else
                return; 
            end
        end
    end

end