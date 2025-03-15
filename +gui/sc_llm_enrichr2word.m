function sc_llm_enrichr2word(selpath, parentfig)

    if nargin<2, parentfig = []; end

    if nargin < 1
        selpath = uigetdir;
    end     
    if isempty(selpath), return; end
    if ~isfolder(selpath), return; end

    files = dir(fullfile(selpath, '*_DE_*.xlsx'));
    fileNames1 = string({files(~[files.isdir]).name});
    files = dir(fullfile(selpath, '*_DV_*.xlsx'));
    fileNames2 = string({files(~[files.isdir]).name});
    
    listItems = [fileNames1'; fileNames2']; 

    if gui.i_isuifig(parentfig)
        [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, 'Select Excel Files:');
    else
        [selectedIndex, ok] = listdlg('PromptString', 'Select Excel Files:', ...
                              'SelectionMode', 'multiple', ...
                              'ListString', listItems, ...
                              'ListSize', [260 300], ...
                              'InitialValue', 1:numel(listItems));
    end    
    
    if ok
        selectedfiles = listItems(selectedIndex);
    else
        return;
    end
    
    import mlreportgen.dom.*
    
    fw = gui.myWaitbar;

    for k = 1:length(selectedfiles)
        gui.myWaitbar(parentfig, fw, false, '', ...
            sprintf('%s', selectedfiles(k)), ...
            (k-0.5)/length(selectedfiles));
        infile = fullfile(selpath, selectedfiles(k));
        [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(infile);
        % assignin("base","TbpUp",TbpUp);
        % assignin("base","TmfUp",TmfUp);
        % assignin("base","TbpDn",TbpDn);
        % assignin("base","TmfDn",TmfDn);
        [~, wordfilename] = fileparts(selectedfiles(k));
        [done, outfile] = llm.e_DETableSummary(TbpUp, TmfUp, TbpDn, ...
            TmfDn, wordfilename);


        % files = dir(fullfile(selpath, '*_DP_*.xlsx'));
        % fileNames = string({files(~[files.isdir]).name});
        % doc = Document(outfile, 'docx');
        % open(doc);
        % para = Paragraph("AI generated text");
        % append(doc, para);
        % close(doc);


        if done, rptview(outfile, 'docx'); end
    end
    gui.myWaitbar(parentfig, fw);  
   
    function [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(excelfile)
        TbpUp = [];
        TmfUp = [];
        TbpDn = [];
        TmfDn = [];
        sheetList = sheetnames(excelfile);
        
        sheetToRead = 'Up_250_GO_BP';
        if any(strcmp(sheetList, sheetToRead))
            TbpUp = readtable(excelfile, 'Sheet', sheetToRead);
        end

        sheetToRead = 'Up_250_GO_MF';
        if any(strcmp(sheetList, sheetToRead))
            TmfUp = readtable(excelfile, 'Sheet', sheetToRead);
        end

        sheetToRead = 'Dn_250_GO_BP';
        if any(strcmp(sheetList, sheetToRead))
            TbpDn = readtable(excelfile, 'Sheet', sheetToRead);
        end

        sheetToRead = 'Dn_250_GO_MF';
        if any(strcmp(sheetList, sheetToRead))
            TmfDn = readtable(excelfile, 'Sheet', sheetToRead);
        end        
    end
end