function sc_llm_enrichr2word(selpath)

    if nargin < 1, selpath = uigetdir; end
    
    files = dir(fullfile(selpath, '*_DE_*.xlsx'));
    fileNames1 = string({files(~[files.isdir]).name});
    files = dir(fullfile(selpath, '*_DV_*.xlsx'));
    fileNames2 = string({files(~[files.isdir]).name});
    
    listItems = [fileNames1'; fileNames2']; 
    [selectedIndex, ok] = listdlg('PromptString', 'Select Excel Files:', ...
                                  'SelectionMode', 'multiple', ...
                                  'ListString', listItems, ...
                                  'ListSize', [260 300], ...
                                  'InitialValue', 1:numel(listItems));
    if ok
        selectedfiles = listItems(selectedIndex);
    else
        return;
    end
    
    import mlreportgen.dom.*
    fw = gui.gui_waitbar_adv;              
    for k = 1:length(selectedfiles)
        gui.gui_waitbar_adv(fw, (k-0.5)/length(selectedfiles), ...
            sprintf('%s', selectedfiles(k));
        infile = fullfile(selpath, selectedfiles(k));
        [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(infile);
        % assignin("base","TbpUp",TbpUp);
        % assignin("base","TmfUp",TmfUp);
        % assignin("base","TbpDn",TbpDn);
        % assignin("base","TmfDn",TmfDn);
        [~, wordfilename] = fileparts(selectedfiles(k));
        [done, outfile] = gui.e_llmsummarizer(TbpUp, TmfUp, TbpDn, TmfDn, wordfilename);
        if done, rptview(outfile, 'docx'); end
    end
    gui.gui_waitbar_adv(fw);    
   
    function [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(excelfile)
        TbpUp = readtable(excelfile, 'Sheet', 'Up_250_GO_BP');
        TmfUp = readtable(excelfile, 'Sheet', 'Up_250_GO_MF');
        TbpDn = readtable(excelfile, 'Sheet', 'Dn_250_GO_BP');
        TmfDn = readtable(excelfile, 'Sheet', 'Dn_250_GO_MF');
    end
end