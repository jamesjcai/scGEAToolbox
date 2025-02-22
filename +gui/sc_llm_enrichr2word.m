function sc_llm_enrichr2word(selpath)

% cd 'C:\Users\jcai\Downloads\scgeatool_DEVPAnalysis_Batch_workingfolder'
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
fw = gui.gui_w              
for k = 1:2  % length(selectedfiles)
    infile = fullfile(selpath, selectedfiles(k));
    [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(infile);

    % assignin("base","TbpUp",TbpUp);
    % assignin("base","TmfUp",TmfUp);
    % assignin("base","TbpDn",TbpDn);
    % assignin("base","TmfDn",TmfDn);

    [~, wordfilename] = fileparts(selectedfiles(k));
    [done, outfile] = gui.e_llmsummarizer(TbpUp, TmfUp, TbpDn, TmfDn, wordfilename);
    if done
        rptview(outfile, 'docx');
    end
end

end


function [TbpUp, TmfUp, TbpDn, TmfDn] = in_gettables(excelfile)
    TbpUp = in_readexceltable(excelfile, 'Up_250_GO_BP');
    TmfUp = in_readexceltable(excelfile, 'Up_250_GO_MF');
    TbpDn = in_readexceltable(excelfile, 'Dn_250_GO_BP');
    TmfDn = in_readexceltable(excelfile, 'Dn_250_GO_MF');
end

function [T] = in_readexceltable(excelfile, sheetname)
    T = readtable(excelfile,'Sheet', sheetname);
    %T = t(:, [3 7]);
    %T.Properties.VariableNames={'Function Term','Genes'};
    %T.Genes = strrep(T.Genes,',', ', ');
end
