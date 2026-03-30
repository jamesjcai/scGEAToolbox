function sc_llm_dp2word(selpath, parentfig)
% SC_LLM_DP2WORD - Generate LLM Word reports from DP batch analysis Excel files
%
% Scans a directory for *_DP_*.xlsx files produced by callback_DPGene2GroupsBatch,
% reads the Up-regulated and Down-regulated program sheets, and calls
% llm.e_DPTableSummary to generate a Word document summary for each file.
%
% Unlike DE/DV reports, no Enrichr step is needed — DP results are already
% at the pathway/program level and can be passed directly to the LLM.
%
% See also: sc_llm_enrichr2word, llm.e_DPTableSummary

if nargin < 2, parentfig = []; end
if nargin < 1
    selpath = uigetdir;
end
if isempty(selpath) || selpath == 0, return; end
if ~isfolder(selpath), return; end

files = dir(fullfile(selpath, '*_DP_*.xlsx'));
listItems = string({files(~[files.isdir]).name});

if isempty(listItems)
    fprintf('No *_DP_*.xlsx files found in %s\n', selpath);
    return;
end

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, ...
        'Select Excel Files:', listItems);
else
    [selectedIndex, ok] = listdlg('PromptString', 'Select Excel Files:', ...
        'SelectionMode', 'multiple', ...
        'ListString', listItems, ...
        'ListSize', [260 300], ...
        'InitialValue', 1:numel(listItems));
end

if ~ok, return; end
selectedfiles = listItems(selectedIndex);

fw = gui.myWaitbar(parentfig);

for k = 1:length(selectedfiles)
    gui.myWaitbar(parentfig, fw, false, '', ...
        sprintf('%s', selectedfiles(k)), ...
        (k-0.5)/length(selectedfiles));

    infile = fullfile(selpath, selectedfiles(k));
    [~, wordfilename] = fileparts(selectedfiles(k));

    % Read Up-regulated and Down-regulated sheets written by sc_dpg batch
    Tup = [];
    Tdn = [];
    try
        sheetList = sheetnames(infile);
        if any(strcmp(sheetList, 'Up-regulated'))
            Tup = readtable(infile, 'Sheet', 'Up-regulated');
        end
        if any(strcmp(sheetList, 'Down-regulated'))
            Tdn = readtable(infile, 'Sheet', 'Down-regulated');
        end
    catch ME
        fprintf('Could not read %s: %s\n', selectedfiles(k), ME.message);
        continue;
    end

    [done, outfile] = llm.e_DPTableSummary(Tup, Tdn, wordfilename, selpath);
    if done, rptview(outfile, 'docx'); end
end

gui.myWaitbar(parentfig, fw);
end
